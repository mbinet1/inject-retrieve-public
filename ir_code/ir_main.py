import all_methods
import numpy as np 
import current_parameters as prm
import time

ir = all_methods.InjectRetrieve()

ir.BeforeLoop()
while ir.count_samples < prm.nb_samples:
    ir.BeforeModel()
    if ir.CarryLoop():
        ir.CreateModel()
        ir.Injection()
        ir.FlatteningRemoval()
        ir.TransitSearch()
        ir.MeanDepth()
        if prm.is_notebook:
            ir.ShowCurves([0,1,2])
            ir.FoldedCurve(True)
        ir.ComputeResults()
        ir.Outcomes()
        ir.WriteResults()
        ir.LoopOutputs()
ir.AfterLoop()