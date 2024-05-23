#!/usr/bin/env nextflow

params.greeting = 'Nextflow basics: testing parameters'
input_ch = Channel.of(params.greeting)

process SPLITLETTERS {
    input: 
    val x

    /* A process must have at most one input block, 
    inputs can be: 
    val: values,
    path: path to files and directories, paths are usually reolved relative to the execution path.
    env: Set an environment variable in the process script
    stdin: Forward the input from stdin to the process
    tuple: Handle a group of input values of any qualifiers.
    each: Execute the process for each element in the input parameters.
    */

    output:
    path 'chunk_*' // Tells the process to expect an output file(s) path with a filename starting with 'chunk_'. The process sends the output as a channel.

    script: 
    /* The script to execute; three double quotes define a code block.*/
    """
    printf "$x" | split -b 6 - chunk_
    """ 
}

process CONVERTTOUPPER {
    input:
    path y //Assign the input for this process to the 'y' variable (path type qualifier).

    output:
    stdout // Tells the process to expect output as standard output (stdout) and send it as a channel. 

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}


workflow { // Start the workflow scope where each process can be called.
   letters_ch = SPLITLETTERS(input_ch) // Execute the SPLITLETTERS process on the 'input_ch' input channel and store the results into the 'letters_ch' channel
   results_ch = CONVERTTOUPPER(letters_ch.flatten()) // Convert to upper case and flatten the output using the 'flatten()' operator which transforms the input channel such that every item is a separate element. Store the output into the 'results_ch' channel
   results_ch.view { it } // Print the results using the 'view' operator
}

// Next: execute the script from the command-line: nextflow run basic.nf