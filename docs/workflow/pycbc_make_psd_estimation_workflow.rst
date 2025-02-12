=========================================================
Estimating average noise PSDs over long stretches of data
=========================================================

It can be useful to estimate the average noise PSD of a long period of data, for
instance for building template banks manually or doing bank simulations. The
program ``pycbc_make_psd_estimation_workflow`` is devoted to this task.

The program sets up a Pegasus workflow which is basically a small subset of a
coincident search workflow:

 * Find the data frames and the segments
 * Segment the analyzable data in each detector
 * Run ``pycbc_calculate_psd`` to estimate the PSD in each segment in each detector
 * Run ``pycbc_average_psd`` to combine the PSD estimates over time for each detector, as well as over time and detectors.

Configuration file
------------------

``pycbc_make_psd_estimation_workflow`` is configured through an .ini file,
similarly to search workflows. An example for ER8 data, broken by sections, is
given here.

::

    [workflow]
    start-time = 1123858817
    end-time = 1125217722
    h1-channel-name = H1:GDS-CALIB_STRAIN
    l1-channel-name = L1:GDS-CALIB_STRAIN
    file-retention-level = no_intermediates

    [workflow-ifos]
    h1 =
    l1 =

    [workflow-datafind]
    datafind-h1-frame-type = H1_HOFT_C00
    datafind-l1-frame-type = L1_HOFT_C00
    datafind-method = AT_RUNTIME_SINGLE_FRAMES
    datafind-check-segment-gaps = update_times
    datafind-check-frames-exist = raise_error
    datafind-check-segment-summary = warn

    [workflow-segments]
    segments-h1-science-name = H1:DMT-ANALYSIS_READY:1
    segments-l1-science-name = L1:DMT-ANALYSIS_READY:1
    segments-database-url = https://segments.ligo.org
    segments-veto-definer-url = file:///home/tito/er8/H1L1V1-ER8_CBC_OFFLINE.xml
    segments-science-veto = 1
    segments-veto-groups =
    segments-final-veto-group = 12H
    segments-method = ALL_SINGLE_IFO_TIME

    [datafind]
    urltype = file

    [workflow-matchedfilter]
    matchedfilter-method = WORKFLOW_INDEPENDENT_IFOS
    analysis-length = 2048
    min-analysis-segments = 15
    max-analysis-segments = 15
    output-type = hdf

The above sections control which data are being used and they are basically the
same as found in a coincidence workflow. In this example, ER8 data are used and
the analyzable time is broken into 2048 s segments for PSD estimation.

::

    [executables]
    segment_query = ${which:ligolw_segment_query_dqsegdb}
    segments_from_cats = ${which:ligolw_segments_from_cats_dqsegdb}
    llwadd = ${which:ligolw_add}
    ligolw_combine_segments = ${which:ligolw_combine_segments}
    plot_segments = ${which:pycbc_page_segments}
    calculate_psd = ${which:pycbc_calculate_psd}
    average_psd = ${which:pycbc_average_psd}
    plot_spectrum = ${which:pycbc_plot_psd_file}
    plot_range = ${which:pycbc_plot_range}
    page_segtable = ${which:pycbc_page_segtable}
    page_segplot = ${which:pycbc_page_segplot}

The above section specifies the location of the various executables called by
the workflow. The ``${which:X}`` syntax replaces the line with the full path to
the executable, wherever that happens to be at the time of running
``pycbc_make_psd_estimation_workflow``.

::

    [llwadd]

    [segments_from_cats]

    [ligolw_combine_segments]

I'm not actually sure these sections are needed: FIXME ;)

::

    [calculate_psd]
    cores = 4
    low-frequency-cutoff = 10
    pad-data = 8
    strain-high-pass = 8
    sample-rate = 4096
    segment-length = 256
    segment-start-pad = 64
    segment-end-pad = 64
    psd-estimation = median
    psd-segment-length = 16
    psd-segment-stride = 8

    [calculate_psd-h1]
    channel-name = H1:GDS-CALIB_STRAIN

    [calculate_psd-l1]
    channel-name = L1:GDS-CALIB_STRAIN

    [pegasus_profile-calculate_psd]
    condor|request_cpus = 4

The above sections control how the PSD is estimated in each segment. The program
devoted to this is ``pycbc_calculate_psd``, see its ``--help`` for details. In
this example, two instances of ``pycbc_calculate_psd`` are launched (one per
detector) and each instance uses 4 CPU cores. For details on PSD estimation,
see for instance the `FindChirp paper <http://arxiv.org/abs/gr-qc/0509116>`_.

::

    [average_psd]

The above section controls how the averaging of the PSDs over time and detector
is done, i.e. it contains options for the ``pycbc_average_psd`` program.
Currently the program does not take options and the only supported averaging
method is the harmonic mean.

::

    [plot_segments]

    [plot_range]
    mass1 = 1.4
    mass2 = 1.4
    approximant = SPAtmplt

    [plot_spectrum]
    psd-model = aLIGOZeroDetHighPower

    [page_segtable]

    [page_segplot]

The above sections control plotting jobs.

Generating and running the workflow
-----------------------------------

Once you have an .ini file at ``/path/to/ini/file``, create the workflow in the
following way:

::

    pycbc_make_psd_estimation_workflow \
        --workflow-name RUN_NAME \
        --output-dir /path/to/run/directory \
        --config-files /path/to/ini/file

``RUN_NAME`` should be replaced with a meaningful descriptive name for the
workflow and ``/path/to/run/directory`` should point to the directory where the
run is supposed to take place. Once the workflow is generated, move to
``/path/to/run/directory`` and start the workflow with

::

    pycbc_submit_dax \
        --dax RUN_NAME.dax \
        --accounting-group ACCOUNTING_TAG

where again ``RUN_NAME`` and ``ACCOUNTING_TAG`` should be given meaningful
values. When the workflow completes, the average PSDs should be available in
``/path/to/run/directory/psds`` and diagnostic plots should be in
``/path/to/run/directory/plots``.
