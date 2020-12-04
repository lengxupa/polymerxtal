#!/usr/bin/env python

import os

siteDesignator = "newton"
monitorRoot = os.path.join(os.sep, 'home', 'HUBzero', 'submit')
qstatCommand = "qstat -u HUBzero"
monitorLogLocation = "logs"
monitorLogFileNme = "monitorPBS.log"
historyFileName = "monitorPBS.history"

SITEDESIGNATOR = "clusterPBS"
MONITORROOT = os.path.join(os.sep, 'home', 'yourhub', 'submit')
QSTATCOMMAND = "/usr/pb/bin/qstat -u yourhub"
MONITORLOGLOCATION = os.path.join(MONITORROOT,'logs')

