import os, sys


def fchdir(fd):
    # initial directory
    cwd = os.getcwd()

    # trying to insert to fd directory
    try:
        os.chdir(fd)
        # print("Inserting inside-", os.getcwd())

    # Caching the exception
    except:
        print("Something wrong with specified directory. Exception- ", sys.exc_info())
        return -1

        # handling with finally
        # finally:
        print("Restoring the path")
        os.chdir(cwd)
        print("Current directory is-", os.getcwd())
