# ============================================================================
# rng.py -- Mersenne Twister pseudorandom number generator
# ----------------------------------------------------------------------------
# Author (wrapper): Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================
# M. Matsumoto and T. Nishimura, "Mersenne Twister: a 623-dimensionally
# equidistributed uniform pseudo-random number generator", ACM Transactions
# on Modeling and Computer Simulation, Vol 8, No. 1, Jan 1998.
# doi:10.1145/272991.272995
# ----------------------------------------------------------------------------
# The following implementation of the Mersenne Twister random number generator
# is adapted from the 2002 version released by the authors, Nishimura and
# Matsumoto.  The copyright notice distributed with their code follows:
#
# Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#
#   3. The names of its contributors may not be used to endorse or promote
#      products derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Any feedback is very welcome.
# http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
# email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
# ============================================================================

from .utils import FREE

# Mersenne Twister definitions
MT_M = 397
MT_N = 624
MT_UPPER_MASK = 0x80000000  # most significant w-r bits
MT_LOWER_MASK = 0x7fffffff  # least significant r bits
MT_MATRIX_A = 0x9908b0df  # constant vector a


# State
class RNG:
    def __init__(self):
        self.mt = []
        self.mag01 = []
        self.mti = 0
        self.use_normal_value = 0
        self.normal_value = 0.
