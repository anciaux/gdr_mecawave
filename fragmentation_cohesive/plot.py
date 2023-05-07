

def addPlotVtk(model):
    model.setBaseName("bar")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("velocity")
    model.addDumpFieldVector("internal_force")
    model.addDumpFieldVector("acceleration")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")
    model.addDumpField("material_index")

    # VTK plot for Cohesive model
    model.setBaseNameToDumper("cohesive elements", "cohesive")
    model.addDumpFieldVectorToDumper(
        "cohesive elements", "displacement")
    model.addDumpFieldToDumper("cohesive elements", "damage")
    model.addDumpFieldVectorToDumper(
        "cohesive elements", "tractions")
    model.addDumpFieldVectorToDumper("cohesive elements", "opening")


def addVtkFiles(model, n_steps_to_save):
    model.dump()
    model.dump("cohesive elements")
