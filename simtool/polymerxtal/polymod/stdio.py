class FILE:
    def __init__(self):
        self.path = ""

    def printf(self, format):
        if self.path:
            self.file.write(format)
        else:
            print(format)

    def fclose(self):
        if self.path:
            self.file.close()

    def rewind(self):
        self.file.seek(0)

    def readline(self):
        return self.file.readline()


stderr = FILE()
stdout = FILE()


def fputc(character, stream):
    stream.printf(character)


def fgets(length, stream):
    return stream.file.readline()
