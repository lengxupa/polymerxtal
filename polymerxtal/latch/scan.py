from .lex import yylex, yyget_text, yyget_lineno


class Scanner:
    def __init__(self, path):
        self.path = path
        self.f = open(path, "r")
        self.scanner = self.f.readline()
        self.lineno = 1
        self.toklen = 0
        self.reuse = 0

    def getToken(self):
        if self.reuse:
            self.reuse = 0
        else:
            self.tokval = yylex(self.scanner)
            self.tokstr = yyget_text(self.scanner)
            self.toklen = yyget_leng(self.scanner);
