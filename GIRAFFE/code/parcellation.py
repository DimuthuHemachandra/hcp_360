from nipype.interfaces.base import (
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    File
)
import os
#import traits

class hcp_360_InputSpec(CommandLineInputSpec):
    script = File(exists=True, argstr='%s', mandatory=True, position=0, desc='bash_file')
    input_file = File(argstr='-L %s -a HCP-MMP1 -f 1 -l 1', mandatory=True, desc='subject_file') #You can have optional flags here at the end too.
    #input_file = File(desc="File", exists=True, mandatory=True, argstr="-L %s -a HCP-MMP1")
    SUBJECTS_DIR = File(argstr='-k %s', mandatory=True, desc='Freesurfer $SUBJECTS_DIR')
    OUT_DIR = File(argstr='-d %s', mandatory=True, desc='Output directory')
    #parameter = traits.Int(desc = "some parameter", argstr = "--k %d")



class hcp_360_OutputSpec(TraitedSpec):
    output_file = File(desc = "out file", exists = True)

class hcp_360(CommandLine):
    input_spec = hcp_360_InputSpec
    output_spec = hcp_360_OutputSpec
    #_cmd = 'bash /home/ylu/Desktop/hcp_360/GIRAFFE/code/derivatives/freesurfer/create_subj_volume_parcellation.sh'
    _cmd = 'bash'


"""
if __name__ == '__main__':

    #os.system('export SUBJECTS_DIR=~/Desktop/sub-CT01/derivatives/freesurfer15/freesurfer/')
    zipper = GZipTask(input_file='subject.txt', SUBJECTS_DIR='/home/ylu/Desktop/sub-CT01/derivatives/freesurfer15/freesurfer/', OUT_DIR='output/')
    print(zipper.cmdline)
    zipper.run()
    #print(zipper._list_outputs())"""
