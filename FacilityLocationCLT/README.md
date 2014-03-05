Facility Location CLT
====

This example generates a small, snythetic facility location example and compares three data-driven methods for solving it:
 - Fixing uncertain demand at its nominal value
 - Naively fitting a CLT style uncertainty set
 - Tuning a normalized CLT uncertainty set
 
The file facilityLocation.jl does most of the work of simulation, and solving using Julia/JuMP.
The file genPlots.r is a simple R script to create some plots.  
