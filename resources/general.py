from inspect import getsourcefile
from os.path import abspath

def locate_dir(folder_name, root=None):
    # does not handle length error...
    curr_dir = abspath(getsourcefile(lambda:0))

    if root and type(root) == int:
        split = curr_dir.split("\\")[0:-root]
        return '\\'.join(split) + f'\\{folder_name}'

    elif root == None:
        return curr_dir + f'\\{folder_name}'

    else:
        raise TypeError('Root location must be int.')