# ============================================================================
# timer.py -- Timer functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import time


class Timer:
    def __init__(self):
        self.t0 = 0.0
        self.tv = 0.0

    # ============================================================================
    # startTimer()
    # ----------------------------------------------------------------------------
    # Result: mark a starting time
    # ============================================================================
    def startTimer(self):
        self.t0 = time.time()

    # ============================================================================
    # getSeconds()
    # ----------------------------------------------------------------------------
    # Result: return the current time to microseconds
    # ============================================================================
    def getSeconds(self):
        self.tv = time.time()
        return self.tv

    # ============================================================================
    # getElapsedTime()
    # ----------------------------------------------------------------------------
    # Result: return the time in seconds, to microseconds, since startTimer()
    # call
    # ============================================================================
    def getElapsedTime(self):
        return self.getSeconds() - self.t0
