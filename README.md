# Winograd
Matrix generation applications for Winograd family of convolution algorithms
In the project there are four applications:

1. winograd_modified
   It creates transformation matrices for Winograd convolution.
   Usage:
     bin/winograd_modified <kernel_size> <output_size> N x <point>
     where:
       number of points N = kernel_size + output_size - 2
       and point is pair of integers: (nominator, denominator)
   In creates three transformation matrices:
     - A^T and saves it in: matrix_output/ATMat.txt
     - B^T and saves it in: matrix_output/BTMat.txt
     - G and saves it in: matrix_output/GMAT.txt
   Example:
     bin/winograd_modified 3 4 0 1 -1 2 -1 1 1 1 2 1
     for:
       kernel size = 3
       output size = 4
       and 5 points: 0, -1/2, -1, 1, 2

2. conv_Toom-Cook_1D
   In runs number of simulations using transformation matrices generated
   by winograd_modified. It assumes transformation matrices for matching
   kernel and output size are located in matrix_input directory and named
   GMat.txt, ATMat.txt and BTMat.txt.

   Usage:
     bin/conv_Toom-Cook_1D <kernel_size> <output_size> [no_experiments]
     where:
       kernel_size, is size of kernel
       output_size, is size of output
       no_experiments, is number of experiments to run, default is 5000

   Example:
     bin/conv_Toom-Cook_1D 3 4 5

   Example output:
     Experiment parameters:
       Kernel size: 3
       Output size: 4
       Input size:  6
     Performing custom number of experiments: 5
     Experiment results:
       Winograd_huffman:     0.000000193693
       Mixed + Huffman:      0.000000150327
       Winograd:             0.000000176936
       Direct convolution:   0.000000055403

3. conv_Toom-Cook_2D_Huffman
   It runs fixed number of simulations using transformation matrices generated
   by winograd_modified. It assumes transformation matrices for matching
   kernel and output size are located in matrix_input directory and named
   GMat.txt, ATMat.txt and BTMat.txt.

   Usage:
     bin/conv_Toom-Cook_2D_Huffman <kernel_size> <output_size> [no_experiments]
     where:
       kernel_size, is size of kernel
       output_size, is size of output
       no_experiments, is number of experiments to run, default is 5000

   Example:
     bin/conv_Toom-Cook_2D_Huffman 3 4 5

   Example output:
     Experiment parameters:
       Kernel size: 3
       Output size: 4
       Input size:  6
     Performing custom number of experiments: 5
     Experiment results:
       Winograd mixed Huffman:       0.000004165826
       Winograd:                     0.000004224974
       Winograd huffman:             0.000004022474
       Direct:                       0.000000783162

4. winograd_modified_full
   It creates full Winograd transformation matrices based on transformation
   matrices. Transformation matrices are located in matrix_input directory.
   New transforamtion matrices are out in matrix_output directory.

   Usage:
     bin/winograd_modified_full

   Example:
     bin/winograd_modified_full
     It will take GMat_2.txt, ATMat_2.txt and BTMat_2.txt from matrix_input
     and will create GMat_Win_6_sym.txt, ATMat_Win_6_sym.txt and
     BTMat_Win_6_sym.txt files with transformation matrices.
