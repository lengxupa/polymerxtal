class FILE:
    def __init__(self):
        self._ptr = ''
        self._cnt = 0
        self._base = ''
        self._flag = 0
        self._file = 0
        self._charbuf = 0
        self._bufsiz = 0
        self._tmpfname = ''


def rewind(f):
    f.seek(0)
