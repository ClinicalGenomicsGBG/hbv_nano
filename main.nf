#!/usr/bin/env nextflow

params.greeting = "hbv_nano"

greeting_ch = Channel.of(params.greeting)

workflow {}