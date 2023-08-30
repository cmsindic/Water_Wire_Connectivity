import os

for dir in os.scandir('../'):
    if 'ases' in dir.name[-4:]:
        for obj in os.scandir(dir):
            if os.path.isdir(obj):
                for file in os.scandir(obj):
                    if 'index' in file.name:
                        p = file.path
                        pl = p.split('/')
                        new = os.path.join(pl[1],pl[-1])
                        os.system('cp {} {}'.format(p,new))
