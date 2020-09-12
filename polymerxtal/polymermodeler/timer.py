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


class timeval:
    def __init__(self):
        self.tv_sec = 0.
        self.tv_usec = 0.


class Timer:
    def __init__(self):
        self.t0 = 0.
        self.tv = timeval()

        # ============================================================================
        # startTimer()
        # ----------------------------------------------------------------------------
        # Result: mark a starting time
        # ============================================================================
        def startTimer(self):
            self.t0 = time.time()
