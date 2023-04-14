#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

process pythonTask {
    """
    #!/usr/bin/python3

    x = 'Hello'
    y = 'world!'
    print "%s - %s" % (x,y)
    """
}

workflow {
    pythonTask()
}