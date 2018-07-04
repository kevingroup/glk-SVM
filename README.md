## glk-SVM 

### Introduction

 The glk-SVM provides implementation of a novel SVM kernel method using glk-patterns as features for DNA sequences and the related classification tasks.

### Installation

Download the zip file and  extract the source codes.

`$ cd glk-SVM`

`$ make`

If successful, you could find the following executables in the current directory:

- glk_kernel
- glk_train
- glk_classify

`make install` will simply copy these three executables to the `../bin` direcory

### Usage

Here we introduce the basic workflow of `glk-SVM`.  Please refer to help messages for more detailed information of each program.  You can find the instructions by running the programs without any argument.

####(1) Calculating kernel matrix

A pre-computed kernel matrix is needed to train SVM classifiers. Here we are going to use `test/test_positive.fa` as a positive sequence set, and `test/test_negative.fa` as a negative set. The resulting kernel matrix is stored in `test_kernel.out`.

Type:

`$ ./glk_kernel -g 2 -l 10 -k 8 test/test_positives.fa test/test_negatives.fa test_kernel.out`

The kernel matrix is in `test_kernel.out`. We also add the same file in the test directory `test/test_kernel.out` and these two files should be identical.

#### (2) Training SVM

Now we can train the SVM classifier using the kernel matrix from previous step. This training method will take four arguments: kernel file, positive sequence file, negative sequence file, and output file prefix.

Type:

`./glk_train test_kernel.out test/test_positives.fa test/test_negatives.fa test`

The training process will generate two files, `test_svalpha.out` and `test_svseq.fa`. `test_svalpha.out` contains the weight values for the support vector sequences, `test_svseq.fa` contains the sequences in fasta format. These two files will be used in the next classification step.

#### (3) Classification using SVM

Now we can use `glk_classify` to score/classify any set of DNA sequences. Here, we will use trained classifier to score the input positive training sequences. We should specify the same set of parameters as used in the first step to ensure correct classification. 

Type:

  `$ ./glk_classify -g 2 -l 10 -k 8 test/test_positives.fa test_svseq.fa test_svalpha.out pos_classify.out`

### Distribution

- README    -- This file
- Makefile  -- Makefile for compile
- src/      -- The source code files
- test/     -- Example sequence files for testing



