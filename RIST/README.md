## RIST
RIST is to refine the insertion sites and obtain the target site duplication sequences (TSDs) as well as other informations of retroduplications detected by GRIPper. 
Additionally, the signals extracted by RIST could be applied to confirm the somatic and loss origin of retroduplications for paired tumor-normal samples.

## Inputs
The input of RIST include the file containing the GRIPper output, please see the format in example/demo_gripper_output.txt.
Also the sample information should be provided, see example/demo_sample_info.txt as an example.

In addtion, the directions of samtools, reference genome and bam file etc, should be added to the file config.txt.

#### Running the test:
  cd example
  sh test.sh

Then see the demo results in the `result` folder.

Please feel free to contact sunnyzxh@foxmail.com if you have any questions in running RIST.
