# Copyright (C) 2015  Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module is responsible for setting up plotting jobs.
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""
from pycbc.workflow.core import FileList, makedir, Executable
 
def excludestr(tags, substr):
    if substr is None:
        return tags
    if isinstance(substr, list):
        if len(substr) > 1:
            tags = excludestr(tags, substr[1:])
        substr = substr[0]   
    return [tag for tag in tags if substr not in tag]

def requirestr(tags, substr):
    if substr is None:
        return tags
    return [tag for tag in tags if substr in tag]
 
class PlotExecutable(Executable):
    """ plot executable
    """
    current_retention_level = Executable.FINAL_RESULT

def make_template_plot(workflow, bank_file, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_bank', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--bank-file', bank_file)
    
    if workflow.cp.has_option_tags('workflow-coincidence', 'background-bins', tags=tags):
        bins = workflow.cp.get_opt_tags('workflow-coincidence', 'background-bins', tags=tags)
        node.add_opt('--background-bins', bins)
    
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files[0]     

def make_range_plot(workflow, psd_files, out_dir, exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_range'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_range', ifos=workflow.ifos,
                              out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_list_opt('--psd-files', psd_files)
        node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_spectrum_plot(workflow, psd_files, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_spectrum', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--psd-files', psd_files)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files[0]
 
def make_segments_plot(workflow, seg_files, out_dir, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_segments', ifos=workflow.ifos,
                         out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
        
def make_foreground_table(workflow, trig_file, bank_file, ftag, out_dir, 
                          singles=None, extension='.html', tags=None):
    if tags is None:
        tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_foreground', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--bank-file', bank_file)
    node.add_opt('--foreground-tag', ftag)
    node.add_input_opt('--trigger-file', trig_file)
    if singles is not None:
        node.add_input_list_opt('--single-detector-triggers', singles)
    node.new_output_file_opt(bank_file.segment, extension, '--output-file')
    workflow += node
    return node.output_files[0]

def make_sensitivity_plot(workflow, inj_file, out_dir, exclude=None, require=None, tags=[]):
    makedir(out_dir)   
    secs = requirestr(workflow.cp.get_subsections('plot_sensitivity'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_sensitivity', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_coinc_snrchi_plot(workflow, inj_file, inj_trig, stat_file, trig_file, out_dir,
                           exclude=None, require=None, tags=[]):
    makedir(out_dir)    
    secs = requirestr(workflow.cp.get_subsections('plot_coinc_snrchi'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_coinc_snrchi', ifos=inj_trig.ifo,
                    out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_opt('--found-injection-file', inj_file)
        node.add_input_opt('--single-injection-file', inj_trig)
        node.add_input_opt('--coinc-statistic-file', stat_file)
        node.add_input_opt('--single-trigger-file', trig_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_inj_table(workflow, inj_file, out_dir, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_injections', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--injection-file', inj_file)
    node.new_output_file_opt(inj_file.segment, '.html', '--output-file')
    workflow += node   

def make_seg_table(workflow, seg_files, seg_names, out_dir, tags=None):
    """ Creates a node in the workflow for writing the segment summary
    table. Returns a File instances for the output file.
    """
    seg_files = list(seg_files)
    seg_names = list(seg_names)
    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_segtable', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    node.add_opt('--segment-names', ' '.join(seg_names))
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_veto_table(workflow, out_dir, vetodef_file=None, tags=None):
    """ Creates a node in the workflow for writing the veto_definer
    table. Returns a File instances for the output file.
    """
    if vetodef_file is None:
        vetodef_file = workflow.cp.get_opt_tags("workflow-segments",
                                           "segments-veto-definer-file", [])
    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_vetotable', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_opt('--veto-definer-file', vetodef_file)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_seg_plot(workflow, seg_files, out_dir, seg_names=None, tags=None):
    """ Creates a node in the workflow for plotting science, and veto segments.
    """

    seg_files = list(seg_files)
    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_segplot', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    node.add_opt('--segment-names', ' '.join(seg_names))
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_ifar_plot(workflow, trigger_file, out_dir, tags=None):
    """ Creates a node in the workflow for plotting cumlative histogram
    of IFAR values.
    """

    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_ifar', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--trigger-file', trigger_file)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files[0]

def make_snrchi_plot(workflow, trig_files, veto_file, veto_name, 
                     out_dir, exclude=None, require=None, tags=[]):
    makedir(out_dir)    
    secs = requirestr(workflow.cp.get_subsections('plot_snrchi'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        for trig_file in trig_files:
            node = PlotExecutable(workflow.cp, 'plot_snrchi',
                        ifos=trig_file.ifo, 
                        out_dir=out_dir, 
                        tags=[tag] + tags).create_node()

            node.set_memory(15000)
            node.add_opt('--segment-name', veto_name)
            node.add_input_opt('--trigger-file', trig_file)
            node.add_input_opt('--veto-file', veto_file)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node  
            files += node.output_files
    return files

def make_foundmissed_plot(workflow, inj_file, out_dir, exclude=None, require=None, tags=[]):
    makedir(out_dir)   
    secs = requirestr(workflow.cp.get_subsections('plot_foundmissed'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        exe = PlotExecutable(workflow.cp, 'plot_foundmissed', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags)
        node = exe.create_node()        
        ext = '.html' if exe.has_opt('dynamic') else '.png'
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, ext, '--output-file')
        workflow += node   
        files += node.output_files
    return files
    
def make_snrifar_plot(workflow, bg_file, out_dir, closed_box=False, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_snrifar', ifos=workflow.ifos,
                out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--trigger-file', bg_file)
    
    if closed_box:
        node.add_opt('--closed-box')
    
    node.new_output_file_opt(bg_file.segment, '.png', '--output-file')
    workflow += node
    return node.output_files[0]
    
def make_results_web_page(workflow, results_dir):
    template_path = 'templates/orange.html'

    out_dir = workflow.cp.get('results_page', 'output-path')
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'results_page', ifos=workflow.ifos,
                out_dir=out_dir).create_node()
    node.add_opt('--plots-dir', results_dir)
    node.add_opt('--template-file', template_path)
    workflow += node

def make_single_hist(workflow, trig_file, veto_file, veto_name, 
                     out_dir, exclude=None, require=None, tags=[]):
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_hist'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_hist',
                    ifos=trig_file.ifo, 
                    out_dir=out_dir, 
                    tags=[tag] + tags).create_node()
                    
        node.add_opt('--segment-name', veto_name)
        node.add_input_opt('--veto-file', veto_file)
        node.add_input_opt('--trigger-file', trig_file)
        node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_singles_plot(workflow, trig_files, bank_file, veto_file, veto_name, 
                     out_dir, exclude=None, require=None, tags=[]):
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_singles'), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        for trig_file in trig_files:
            node = PlotExecutable(workflow.cp, 'plot_singles',
                        ifos=trig_file.ifo, 
                        out_dir=out_dir, 
                        tags=[tag] + tags).create_node()

            node.set_memory(15000)
            node.add_opt('--segment-name', veto_name)
            node.add_input_opt('--bank-file', bank_file)
            node.add_input_opt('--veto-file', veto_file)
            node.add_opt('--detector', trig_file.ifo)
            node.add_input_opt('--single-trig-file', trig_file)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node
            files += node.output_files
    return files
