"""
Line magic to allow the 'use' command to load hub environment modules.
This works on the current environment, so we cannot simply fork a shell.
Maybe some trick to fork a shell and return the modified environment would work?
For now lets see how this works.  Just implement "setenv", "prepend", and "use".
Works for all currently installed modules except two old ones that use local shell variables.
"""

# from IPython.core.magic import register_line_magic
import os
from string import Template
import sys

EPATH = []
d = {}


def check_nanohub():
    if os.environ["HOME"].split("/")[2] == "nanohub":
        global EPATH
        EPATH = os.environ["ENVIRON_CONFIG_DIRS"].split()
        use("lammps-31Mar17")
        use("ovito-2.9.0")
        return True
    else:
        return False


def setenv(line):
    name = line[0]
    val = " ".join(line[1:])
    os.environ[name] = val


def prepend(line):
    global d
    name, val = line
    val = Template(val).safe_substitute(d)
    try:
        oldval = os.environ[name]
        val = "%s:%s" % (val, oldval)
    except:  # noqa: E722
        pass
    os.environ[name] = val

    if name == "PYTHONPATH":
        for p in reversed(val.split(":")):
            sys.path.insert(1, p)


def _set(a, b):
    global d
    t = Template(b)
    b = t.safe_substitute(d)
    d[a] = b


def _use(name):

    fname = None
    for e in EPATH:
        ename = os.path.join(e, name)
        if os.path.isfile(ename):
            fname = ename
            break
    if fname is None:
        raise ValueError("use: could not find module: %s" % name)

    with open(fname) as fp:
        for line in fp:
            sline = line.strip().split()
            if sline == []:
                continue
            if sline[0] == "prepend":
                prepend(sline[1:])
                continue
            if sline[0] == "setenv":
                setenv(sline[1:])
                continue
            if sline[0] == "use":
                _use(sline[-1])
                continue
            line = line.split("=")
            if len(line) == 2:
                _set(line[0].strip(), line[1].strip())


def use(name):
    _use(name)
