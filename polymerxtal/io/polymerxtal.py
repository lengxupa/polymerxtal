"""
Functions for manipulating PolymerXtal files.
"""

import os

from .polymod import generate_polymod_input


def read_input(infile):
    ifile = open(infile)
    args = {}
    for line in ifile.readlines():
        valid_line = line.split("#")[0]
        ln = valid_line.split()
        if len(ln) < 1:
            continue
        args[ln[0]] = ln[1:]
    ifile.close()
    return args


def check_args(args):

    # Check element input
    element_flag = 0
    element_list = []
    if "element" in args:
        if len(args["element"]) > 0:
            num_element = int(eval(args["element"][0]))
            if len(args["element"]) == (1 + num_element):
                element_flag = 1
                element_list = args["element"][1:]
                print(str(num_element) + " elements:")
                print(element_list)
        else:
            print(
                "Please specify number of element types followed by a sequence of the element(s)"
            )
            return False
    if not element_flag:
        return False


def complete_args(args):

    if not check_args(args):
        return False

    # Add output data file name
    output = "polymer_crystal"
    if "output" in args:
        if len(args["output"]) == 1:
            output = args["output"][0]
            print("Output file %s" % output)
        else:
            print("Please define an output name")
            return False
    else:
        print("Default output file %s" % output)

    # Add polymer type information
    polymer_type = "PS-atactic"
    flag_polymer_custom = 0
    polymer_types = [
        "PS-atactic",
        "PS-isotactic",
        "PS-syndiotactic",
        "PE",
        "PMMA-atactic",
        "PMMA-isotactic",
        "PMMA-syndiotactic",
        "PP-atactic",
        "PP-isotactic",
        "PP-syndiotactic",
        "POM",
        "PTFE",
        "PVC",
    ]
    if "polymer_type" in args:
        if len(args["polymer_type"]) == 1:
            polymer_type = args["polymer_type"][0]
            if (
                len(polymer_type.split("-")) == 1
                and polymer_type != "PE"
                and polymer_type != "POM"
                and polymer_type != "PTFE"
                and polymer_type != "PVC"
            ):
                polymer_type += "-atactic"
                print(
                    "Default tacticity: atactic, if you need to specify tacticity, please use syntax: polymer_type {polymer}-isotactic, for example, polymer_type PS-isotactic"
                )
            if polymer_type not in polymer_types:
                if polymer_type == "custom":
                    print(
                        "Please specify an input file of custom polymer you will to build"
                    )
                    return False
                else:
                    print(
                        "Polymer type does not exists in our system, Please check README for all incoporated grain types"
                    )
                    return False
            else:
                print("Polymer type: " + polymer_type)
        elif len(args["polymer_type"]) > 1 and args["polymer_type"][0] == "custom":
            polymer_type = args["polymer_type"][0]
            flag_polymer_custom = 1
            polymer_custom_input = args["polymer_type"][1]
        else:
            print("Please specify a polymer type")
            return False
    else:
        print("Default polymer type: PS-atactic")

    molecular_weight = 0
    if flag_polymer_custom:
        polymer_type_custom = generate_polymod_input(polymer_custom_input)
        if polymer_type_custom:
            os.system("latch/latch in.txt > out.txt")
            polymer_type_custom = read_latchout("out.txt", Dir=polymer_type_custom)
            molecular_weight = writepolymodfile(
                polymer_type_custom, "polymer_types/custom/file.txt"
            )
        else:
            return False
