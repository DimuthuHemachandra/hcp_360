#This is a Nipype generator. Warning, here be dragons.
#!/usr/bin/env python

import sys
import nipype
import nipype.pipeline as pe
import argparse

import nipype.interfaces.io as io
import nipype.interfaces.freesurfer as freesurfer

parser = argparse.ArgumentParser(description='Example skullstripping BIDS App entrypoint script.')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
                    'should be stored. If you are running group level analysis '
                    'this folder should be prepopulated with the results of the'
                    'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                    'Multiple participant level analyses can be run independently '
                    '(in parallel) using the same output_dir.',
                    choices=['participant', 'group'])
parser.add_argument('--participant_label', help='The label(s) of the participant(s) '
                   'that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")


args = parser.parse_args()
bids_dir = args.bids_dir
out_dir = args.output_dir

# Extract sub_ids for SelectFiles
if args.participant_label:
    sub_ids = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    sub_ids = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

#Flexibly collect data from disk to feed into workflows.
io_select_files = pe.Node(io.SelectFiles(templates={'T1':'sub-{sub_id}/anat/sub-{sub_id}_T1w.nii.gz'}), name = 'io_select_files')
io_select_files.inputs.base_directory = bids_dir
io_select_files.iterables = [('sub_id', sub_ids)]

#Wraps the executable command ``recon-all``.
freesurfer_recon_all = pe.Node(interface = freesurfer.ReconAll(), name='freesurfer_recon_all')
freesurfer_recon_all.inputs.directive = "all"

#Generic datasink module to store structured outputs
io_data_sink = pe.Node(interface = io.DataSink(), name='io_data_sink')
io_data_sink.inputs.base_directory = out_dir
io_data_sink.inputs.parameterization = True

#Create a workflow to connect all those nodes
analysisflow = nipype.Workflow('MyWorkflow')
analysisflow.connect(io_select_files, "T1", freesurfer_recon_all, "T1_files")
analysisflow.connect(io_select_files, "sub_id", freesurfer_recon_all, "subject_id")
#analysisflow.connect(freesurfer_recon_all, "out_dir", io_data_sink, "recon_results")

#Run the workflow
plugin = 'MultiProc' #adjust your desired plugin here
plugin_args = {'n_procs': 1} #adjust to your number of cores
analysisflow.write_graph(graph2use='flat', format='png', simple_form=False)
analysisflow.run(plugin=plugin, plugin_args=plugin_args)
