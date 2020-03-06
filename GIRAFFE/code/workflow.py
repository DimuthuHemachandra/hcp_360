#!/usr/bin/env python
#This is a Nipype generator. Warning, here be dragons.


import os
import sys
import nipype
import nipype.pipeline as pe
from nipype import Workflow
import argparse
from os.path import join as opj
from shutil import copy

import parcellation as parce

import nipype.interfaces.io as io
import nipype.interfaces.utility as utility

import nipype.interfaces.freesurfer as freesurfer
from nipype.interfaces.utility import IdentityInterface, Function
from nipype import config, logging

config.set('execution', 'remove_unnecessary_outputs', 'false')

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

fs_folder = opj(out_dir, 'freesurfer')  # location of freesurfer folder
os.system('mkdir -p %s'%fs_folder)

#os.system('/neurodocker/startup.sh')

copy("create_subj_volume_parcellation.sh",fs_folder)
copy("lh.HCP-MMP1.annot",fs_folder)
copy("rh.HCP-MMP1.annot",fs_folder)

#Basic interface class generates identity mappings
Parameters = pe.Node(utility.IdentityInterface(fields=["sub_id"]), name='Parameters', iterfield = ['subj_id'])
Parameters.iterables = [('sub_id', sub_ids)]

#Parameters.run()
print(Parameters.outputs)

#Flexibly collect data from disk to feed into workflows.
SelectFiles = pe.Node(io.SelectFiles(templates={'T1':'sub-{sub_id}/anat/sub-{sub_id}_T1w.nii.gz'}), name = 'SelectFiles')
SelectFiles.inputs.base_directory = bids_dir
SelectFiles.inputs.sub_id = "sub_id"
SelectFiles.outputs.sub_id = "sub_id"
#SelectFiles.iterables = [('sub_id', sub_ids)]
SelectFiles.out_dir = out_dir
#NodeHash_1fc3610.iterables = [('sub_id', sub_ids)]

#print(NodeHash_1fc3610.outputs)
############################################################################


#Wraps the executable command ``recon-all``.
freesurfer_recon_all = pe.Node(interface = freesurfer.ReconAll(), name='freesurfer_recon_all')
freesurfer_recon_all.inputs.directive = "all"
freesurfer_recon_all.inputs.subject_id = "sub_ids"
freesurfer_recon_all.inputs.subjects_dir = fs_folder
freesurfer_recon_all.base_dir = opj(out_dir, 'work')

#b = pe.Node(interface=B(), name="b")
#freesurfer_recon_all.iterables = ("subject_id", ["sub-"+s for s in sub_ids])

#Generic datasink module to store structured outputs
io_data_sink = pe.Node(interface = io.DataSink(), name='io_data_sink')
io_data_sink.inputs.base_directory = out_dir
io_data_sink.inputs.subjects_dir = True



################################################################################
def get_paths(subject_id, out_dir, recon_results):
    # Remember that all the necesary imports need to be INSIDE the function for the Function Interface to work!
    import os
    import sys
    import numpy as np
    from shutil import copy

    out_dir = out_dir+"/freesurfer"
    if os.path.exists(out_dir+"/fsaverage"):
    	print("fsaverage already exists")
    else:
    	print("copying fsaverage from freesurfer")
    	os.system('cp $SUBJECTS_DIR/fsaverage '+out_dir)

    if not os.path.exists(out_dir+"/create_subj_volume_parcellation.sh"):
    	print("current wd is ***************",os.getcwd())
    	#copy("/home/ylu/Desktop/sub-CT01/derivatives/freesurfer15/freesurfer/derivatives/work/create_subj_volume_parcellation.sh",out_dir)
      	#copy("/Users/dimuthu/Documents/Robarts/apps/hcp_360/GIRAFFE/code/derivatives/work/create_subj_volume_parcellation.sh",out_dir) 
    
    #copy(out_dir+'/'+subject,out_dir)
    #subject = "sub-"+subject_id
    subject = subject_id

    file_name = out_dir+'/'+subject+'.txt'

    np.savetxt(file_name, [subject], fmt="%s")

    return file_name


#get_paths("CT01", out_dir)

###########################################################################################
hcppaths = pe.Node(Function(function=get_paths, input_names=["subject_id","out_dir","recon_results"],
                            output_names=["txtfilename"]), name="hcppaths")
hcppaths.inputs.out_dir = out_dir
#hcppaths.iterables = ("subject_id", sub_ids)

#hcppaths.run()
###########################################################################################

path_script = fs_folder+'/create_subj_volume_parcellation.sh'

parcel = pe.Node(interface = parce.hcp_360(script=path_script, SUBJECTS_DIR=out_dir+'/freesurfer', OUT_DIR='out_put/'), name='parce')
parcel.inputs.input_file = "input_file"
#parce.inputs.SUBJECTS_DIR = "SUBJECTS_DIR"
#parce.inputs.OUT_DIR = "OUT_DIR"


#parcel.run()
###########################################################################################

#Create a workflow to connect all those nodes
analysisflow = nipype.Workflow('MyWorkflow')
analysisflow.base_dir = out_dir+'/work'
#analysisflow.connect(SelectFiles, "T1", freesurfer_recon_all, "T1_files")

analysisflow.connect(Parameters, "sub_id", SelectFiles, "sub_id")
analysisflow.connect(SelectFiles, "T1", freesurfer_recon_all, "T1_files")
analysisflow.connect(Parameters, "sub_id", freesurfer_recon_all, "subject_id")
analysisflow.connect(Parameters, "sub_id", hcppaths, "subject_id")
analysisflow.connect(freesurfer_recon_all, "T1", hcppaths, "recon_results")
analysisflow.connect(hcppaths, "txtfilename", parcel, "input_file")
#analysisflow.connect(parcel, "output_file", io_data_sink, "base_directory")

#analysisflow.connect(freesurfer_recon_all, "out_dir", io_data_sink, "recon_results")

if __name__ == "__main__":
    #Run the workflow
    plugin = 'MultiProc' #adjust your desired plugin here
    plugin_args = {'n_procs': 8} #adjust to your number of cores
    #analysisflow.write_graph(graph2use='flat', format='png', simple_form=False)
    analysisflow.run(plugin=plugin, plugin_args=plugin_args)
