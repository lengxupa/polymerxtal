# ============================================================================
# os.py -- Operating system specific functions XXX currently Posix only XXX
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import errno, os  # , sys  # noqa: E401

from .unistd import fchdir

# File scope
prevdir = ""


# ============================================================================
# storeDir()
# ----------------------------------------------------------------------------
# Result: store the current working directory
# ============================================================================
def storeDir():
    prevdir = os.getcwd()
    if not prevdir:
        raise IOError("Unable to open current directory: %s" % os.strerror(errno.EIO))


# ============================================================================
# changeDir()
# ----------------------------------------------------------------------------
# Result: change the current working directory; calls choke() on failure
# ============================================================================
def changeDir(directory):
    fd = directory

    if not os.path.exists(fd):
        raise IOError("Unable to open %s: %s" % (directory, os.strerror(errno.ENOENT)))
    if -1 == fchdir(fd):
        raise IOError(
            "Unable to change directory to %s: %s" % (directory, os.strerror(errno.EIO))
        )


# ============================================================================
# restoreDir()
# ----------------------------------------------------------------------------
# Result: change the current working directory back to the value saved by
# storeDir()
# ============================================================================
def restoreDir():
    if (not prevdir) or -1 == fchdir(prevdir):  # noqa: F823
        raise IOError("Unable to return to previous directory")
    prevdir = ""  # noqa: F841
