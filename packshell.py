#!/usr/bin/env python3
'''
Using packmol software to pack molecules.

This programme can pack molecules in a particular structure simply.

It can do these things:
    1. pack molecules as membrane in x, y area in the form of orthorhombic.
    2. pack molecules as membrane in x, y area random.
    3. restrict atoms below or above a line paralelling axles.
    4. More fast than normal packing method.

Packmol: http://www.ime.unicamp.br/~martinez/packmol/

Using option -h can get the usage information.
'''
import argparse
import configparser
import os.path
import time
from collections import OrderedDict
import moltoolkit


def echo_warning(string):
    '''Give a colorful display for string as a warning'''
    warning = '\033[93m'
    endc = '\033[0m'
    print(warning + "Warning: " + string + endc)


def echo_note(string):
    '''Give a colorful display for string as a note'''
    note = '\033[92m'
    endc = '\033[0m'
    print(note + "Note: " + string + endc)


def echo_error(string):
    '''Give a colorful display for string as a fail'''
    fail = '\033[91m'
    endc = '\033[0m'
    print(fail + "Error: " + string + endc)


def trans_num(number):
    '''Translate the number into string.'''
    if isinstance(number, int):
        return str(number) + '.'
    elif isinstance(number, float):
        return str(round(number, 1))
    else:
        return ""


def packmol(inp_file, cmd='packmol', nopause=True):
    '''Pack molecules using packmol software according inp file'''
    exist = False
    cmd = 'packmol'
    for cmdpath in os.environ['PATH'].split(':'):
        if os.path.isdir(cmdpath) and cmd in os.listdir(cmdpath):
            exist = True
    if not exist:
        print("The %s command is not exist in system path. So I sorry not to pack molecules." % cmd)
        exit(1)
    else:
        try:
            if not nopause:
                input("Press Enter to start packing...")
        except KeyboardInterrupt:
            exit(1)
        return os.system('{} < {}'.format(cmd, inp_file))


def ppackmol(inp_file, cores=1, nopause=True):
    '''Pack molecules using packmol software according inp file parallelly'''
    exist = False
    cmd = 'ppackmol'
    for cmdpath in os.environ['PATH'].split(':'):
        if os.path.isdir(cmdpath) and cmd in os.listdir(cmdpath):
            exist = True
    if not exist:
        print("The %s command is not exist in system path. So I sorry not to pack molecules." % cmd)
        exit(1)
    else:
        try:
            if not nopause:
                input("Press Enter to start packing...")
        except KeyboardInterrupt:
            exit(1)
        return os.system('{} {} {}'.format(cmd, cores, inp_file))


class MolPack():
    '''Pack molecules according the given configure file.'''
    def __init__(self, pack_dict):
        '''
        @Args:
            pack_dict: A dict structure including the options. View ReadConfig for details.
        '''
        self.pack = pack_dict
        # Get the 3d size of pdb
        if self.pack['pack_type'] == 'membrane':
            if os.path.exists(self.pack['pdb']):
                # moltoolkit can get the size of the pdb
                mol = moltoolkit.Mol(self.pack['pdb'])
            else:
                echo_error('The file %s is not exist in current directory' % self.pack['pdb'])
                exit(1)
            self.pdb_size = list(mol.lenths_at_axis)
            self.pdb_size.sort()
            echo_note('The 3D size of molecule is {0[0]} {0[1]} {0[2]}. Please check if there is enough space.'.format(self.pdb_size))

    def packmol(self):
        '''Generating the inp file content for Packing molecules.

        @Return
            self.results: a list including complete or partial content of inp file. Because some packing will be divided into several parts, so the former members of self.results are complete content of inp files, which must pack firstly; the last member is the partial content which can be combiled with other packing content.
        '''
        self.results = []
        if self.pack['pack_type'] == 'fixed':
            # fix or move pdb
            self.__fixed_pack()
        elif self.pack['pack_type'] == 'membrane':
            if self.pack['method'] == 'line':
                self.__mem_line_pack()
            elif self.pack['method'] == 'loose':
                self.__mem_loose_pack()
            elif self.pack['method'] == 'copy':
                self.__mem_copy_pack()
            elif self.pack['method'] == 'general':
                self.__mem_general_pack()
            else:
                echo_error("The method must be one of the general, line, loose, copy.")
                exit(1)
        elif self.pack['pack_type'] == 'random':
            self.__random_pack()

        return self.results

    def __fixed_pack(self):
        '''Fix or move this pdb.'''
        self.results.append('structure {0}\n\tnumber 1\n\tfixed {1[0]} {1[1]} {1[2]} {1[3]} {1[4]} {1[5]}\nend structure\n\n'.format(self.pack['pdb'], self.pack['offset']))

    def __mem_general_pack(self):
        '''Packing molecules into a orthorhombic structure. It will restrict position every molecules.
           .       .
               .       .
           .       .
               .       .
           .       .
               .       .
           .       .
               .       .
        '''
        # Remove some molecules
        takeout = self.pack['takeout']
        echo_note('Packing {0[0]}*{0[1]}-{1} molecules in the box of {2[0]} {2[1]} {2[2]} {2[3]} {2[4]} {2[5]}.'.format(self.pack['align_num'], takeout, self.pack['box']))

        box_size = [self.pack['box'][i + 3] - self.pack['box'][i] for i in range(3)]
        pdb_size = self.pdb_size

        self.results.append('')

        y_coord = self.__get_y_coord_iter(self.pack['box'][1], box_size[1], pdb_size[1], self.pack['align_num'][1])
        for (row, y1, y2) in y_coord:
            even = True
            if row % 2 != 0:
                even = False
            x_coord = self.__get_x_coord_iter(self.pack['box'][0], box_size[0], pdb_size[0], self.pack['align_num'][0], even)
            y1 = y1 if y1 - 0.2 < self.pack['box'][1] else y1 - 0.2
            y2 = y2 if y2 + 0.2 > self.pack['box'][4] else y2 + 0.2
            for column, x1, x2 in x_coord:
                if takeout > 0 and column == 0:
                    takeout -= 1
                    continue
                # 稍微讓分子空間寬鬆一點，它才能好好的站整齊
                x1 = x1 if x1 - 0.2 < self.pack['box'][0] else x1 - 0.2
                x2 = x2 if x2 + 0.2 > self.pack['box'][3] else x2 + 0.2
                self.results[0] += 'structure %s\n' % self.pack['pdb']
                # 添加殘基排序选项以使膜堆砌时残基序号递增
                self.results[0] += '\tresnumbers 2 \n'
                self.results[0] += '\tnumber 1\n'
                self.results[0] += '\tinside box'
                z1 = self.pack['box'][2]
                z2 = self.pack['box'][5]
                for m in x1, y1, z1, x2, y2, z2:
                    self.results[0] += ' ' + trans_num(m)
                self.results[0] += '\n'
                self.results[0] += self.__set_panel(x1, y1, z1, x2, y2, z2)
                self.results[0] += 'end structure\n\n'

    def __mem_copy_pack(self):
        '''Packing molecules into a orthorhombic structure. It will first pack one line and then copy this line to other lines.
           .       .
               .       .
           .       .
               .       .
           .       .
               .       .
           .       .
               .       .
        '''
        # Remove some molecules
        echo_note('Packing {0[0]}*{0[1]} molecules in the box of {1[0]} {1[1]} {1[2]} {1[3]} {1[4]} {1[5]}.'.format(self.pack['align_num'], self.pack['box']))

        box_size = [self.pack['box'][i + 3] - self.pack['box'][i] for i in range(3)]
        pdb_size = self.pdb_size

        now = str(time.time()).replace('.', '')
        first_output = 'mid' + now + '.pdb'
        self.results.append('output /tmp/' + first_output + '\n\n')

        y_coord = self.__get_y_coord_iter(self.pack['box'][1], box_size[1], pdb_size[1], self.pack['align_num'][1])

        # Packing the first line
        row, y1, y2 = next(y_coord)
        even = True if row % 2 == 0 else False
        x_coord = self.__get_x_coord_iter(self.pack['box'][0], box_size[0], pdb_size[0], self.pack['align_num'][0], even)
        y1 = y1 if y1 - 0.2 < self.pack['box'][1] else y1 - 0.2
        y2 = y2 if y2 + 0.2 > self.pack['box'][4] else y2 + 0.2
        for column, x1, x2 in x_coord:
            # 稍微讓分子空間寬鬆一點，它才能好好的站整齊
            x1 = x1 if x1 - 0.2 < self.pack['box'][0] else x1 - 0.2
            x2 = x2 if x2 + 0.2 > self.pack['box'][3] else x2 + 0.2
            self.results[0] += 'structure %s\n' % self.pack['pdb']
            # 添加殘基排序选项以使膜堆砌时残基序号递增
            self.results[0] += '\tresnumbers 2 \n'
            self.results[0] += '\tnumber 1\n'
            self.results[0] += '\tinside box'
            z1 = self.pack['box'][2]
            z2 = self.pack['box'][5]
            for m in x1, y1, z1, x2, y2, z2:
                self.results[0] += ' ' + trans_num(m)
            self.results[0] += '\n'
            self.results[0] += self.__set_panel(x1, y1, z1, x2, y2, z2)
            self.results[0] += 'end structure\n\n'

        self.results.append('')
        # The first line
        self.results[1] += 'structure /tmp/%s\n' % first_output
        self.results[1] += '\tnumber 1\n'
        self.results[1] += '\tfixed 0. 0. 0. 0. 0. 0.\n'
        self.results[1] += 'end structure\n\n'
        # Packing other lines
        for row, y1, y2 in y_coord:
            even = True if row % 2 == 0 else False
            self.results[1] += 'structure /tmp/%s\n' % first_output
            self.results[1] += '\tnumber 1\n'
            y1 = y1 if y1 - 0.2 < self.pack['box'][1] else y1 - 0.2
            x_coord = self.__get_x_coord_iter(self.pack['box'][0], box_size[0], pdb_size[0], self.pack['align_num'][0], even)
            column, x1, x2 = next(x_coord)
            x1 = x1 if x1 - 0.2 < self.pack['box'][0] else x1 - 0.2
            self.results[1] += '\tfixed {} {} 0. 0. 0. 0.\n'.format(trans_num(x1), trans_num(y1))
            self.results[1] += 'end structure\n\n'

    def __set_panel(self, *pdb_vol):
        '''Output pannel restriction.

        There is a margin value which is the distance between molecule and pannel. For example, position of a molecule whose x is 10 at most, 5 at least, you want this molecule above x=7, so the margin is 10-7=3; or you want it below x=7, so the margin is 7-5=3.
        '''
        result = ''
        constraint = ('over_x', 'below_x', 'over_y', 'below_y', 'over_z', 'below_z')
        for k in constraint:
            if len(self.pack[k]) == 0:
                continue
            for margin, indices in self.pack[k].items():
                result += '\tatoms'
                for m in indices:
                    result += ' ' + str(m)
                # if the first letter is o, which means over_*
                if k[0] == 'o':
                    result += '\n\t\tover plane '
                    if k[-1] == 'x':
                        result += '1. 0. 0. ' + trans_num(pdb_vol[3] - margin)
                    elif k[-1] == 'y':
                        result += '0. 1. 0. ' + trans_num(pdb_vol[4] - margin)
                    elif k[-1] == 'z':
                        result += '0. 0. 1. ' + trans_num(pdb_vol[5] - margin)
                else:
                    result += '\n\t\tbelow plane '
                    if k[-1] == 'x':
                        result += '1. 0. 0. ' + trans_num(pdb_vol[0] + margin)
                    elif k[-1] == 'y':
                        result += '0. 1. 0. ' + trans_num(pdb_vol[1] + margin)
                    elif k[-1] == 'z':
                        result += '0. 0. 1. ' + trans_num(pdb_vol[2] + margin)
                result += '\n\tend atoms\n'
        return result

    def __get_x_coord_iter(self, init_coord, x_len, pdb_x_len, num, even=True):
        '''Caculate the beginning and end coordinate of molecule at x side.

        @Arg:
            init_coord: the beginning coordinate
            x_len: the length of x side
            pdb_x_len: the length of x side of pdb
            num： numbers of molecules
            even: is this line is a even line

        @Return:
            coord: line number, beginning coordinate, end coordinate
        '''
        if x_len - pdb_x_len * (num + 0.5) < 0:
            echo_warning('The length of x side {0} is not enough to packing {1} molecules whose size is {2}.'.format(x_len, num, pdb_x_len))
        space = (x_len - (num + 0.5) * pdb_x_len) / (num - 0.5)
        if space < 1.5:
            echo_warning('The space %s between two molecules at x side is less than 1.5. It maybe take a long time to pack.' % space)

        coord = init_coord if even else init_coord + 0.5 * (space + pdb_x_len)
        index = 0
        yield (index, coord, coord + pdb_x_len)
        for i in range(1, num):
            coord += pdb_x_len + space
            index += 1
            yield (index, coord, coord + pdb_x_len)

    def __get_y_coord_iter(self, init_coord, y_len, pdb_y_len, num):
        '''Caculate the beginning and end coordinate at y side of a molecule.

        @Arg:
            init_coord: The beginning coordinate at y side.
            y_len: The length of y side of box.
            pdb_y_Len: The length of y side of pdb.
            num：the number of molecules.

        @Return:
            coord: line number, beginning coordinate, end coordinate
        '''
        if y_len - pdb_y_len * num / 2 < 0:
            echo_warning('The length of y side {0} of box is not enough to packing {1} molecules whose size is {2}.'.format(y_len, num, pdb_y_len))

        num_of_line = num // 2 if num % 2 == 0 else num // 2 + 1
        space = (y_len - (num_of_line + 0.5) * pdb_y_len) / (num_of_line - 0.5)
        # If the number of lines at y is not even, add the distance of last line to the space.
        if num % 2 != 0:
            space += (space + pdb_y_len) / (2 * ((num + 1) / 2 - 1))

        # Caculate the distance between two molecules at x side in order to remind user if there is too small space between molecules at random direction.
        x_space = (self.pack['box'][3] - self.pack['box'][0] - (self.pack['align_num'][0] + 0.5) * self.pdb_size[0]) / (self.pack['align_num'][0] - 0.5)
        xy_space = pow(pow(x_space / 2, 2) + pow(space / 2, 2), 0.5)
        if xy_space < 2.5 and x_space < 2.5:
            echo_warning('The space %s between two molecules is less than 2.5. It maybe take a long time to pack.' % xy_space)
        if space < 1.5:
            echo_warning('The space %s between two molecules at y side is less than 1.5. It maybe take a long time to pack.' % space)

        for i in range(0, num):
            # packing y
            # 3   .  init_coord + 0.5 * (space + pdb_y_Len) + space + len
            # 2 .    init_coord + space + pdb_y_len
            # 1   .  init_coord + 0.5 * (space + pdb_len)
            # 0 .    init_coord
            if i % 2 == 0:
                coord = init_coord + i / 2 * (pdb_y_len + space)
            else:
                coord = init_coord + 0.5 * (space + pdb_y_len) + (i - 1) / 2 * (pdb_y_len + space)
            yield (i, coord, coord + pdb_y_len)

    def __mem_line_pack(self):
        '''Packing molecules line by line. I will not restrict every coordinate of molecule but a line of molecules.'''
        echo_note('Packing {0[0]}*{0[1]} molecules in the box of {1[0]} {1[1]} {1[2]} {1[3]} {1[4]} {1[5]}.'.format(self.pack['align_num'], self.pack['box']))

        # box size
        box_size = [self.pack['box'][i + 3] - self.pack['box'][i] for i in range(3)]

        # The 3D size of pdb
        pdb_size = self.pdb_size

        # The following is the algorithm:
        # p is the length of pdb, s is the space of two pdbs, s2 is tail space of the first line or head space of the second line. n is the numbers of one line. x is side-length. The equal is:
        #   p s p s p s p s2 #
        # y .   .   .   .    # n*p + (n-1)s + 0.5*(p+s) = x
        # y   .   .   .   .  #
        # y .   .   .   .    #
        # y   .   .   .   .  #
        # y .   .   .   .    #
        #   x  x   x   x
        # This is a orthorhombic of x*y=4*5.
        self.results.append('')

        y_coord = self.__get_y_coord_iter(self.pack['box'][1], box_size[1], pdb_size[1], self.pack['align_num'][1])
        for (row, y1, y2) in y_coord:
            even = True
            if row % 2 != 0:
                even = False
            x1, x2 = self.__get_x_line_coord(self.pack['box'][0], box_size[0], pdb_size[0], self.pack['align_num'][0], even)
            self.results[0] += 'structure %s\n' % self.pack['pdb']
            # The number of one line molecules
            self.results[0] += '\tnumber ' + str(self.pack['align_num'][0]) + '\n'
            # 添加殘基排序选项以使膜堆砌时残基序号递增
            self.results[0] += '\tresnumbers 2 \n'
            self.results[0] += '\tinside box'
            z1 = self.pack['box'][2]
            z2 = self.pack['box'][5]
            for m in x1, y1, z1, x2, y2, z2:
                self.results[0] += ' ' + trans_num(m)
            self.results[0] += '\n'
            self.results[0] += self.__set_panel(x1, y1, z1, x2, y2, z2)
            self.results[0] += 'end structure\n\n'

    def __get_x_line_coord(self, init_coord, x_len, pdb_x_len, num, even=True):
        '''Caculate the beginning and end coordinate of a line of molecules.

        @Arg:
            init_coord: the beginning coordinate
            x_len: the length of x side
            pdb_x_len: the length of x side of pdb
            num： numbers of molecules
            even: is this line is a even line

        @Return:
            coord: beginning coordinate, end coordinate
        '''
        # Check if there is enough space to packing molecules.
        if x_len - pdb_x_len * (num + 0.5) < 0:
            echo_warning('The length of x side {0} is not enough to packing {1} molecules whose size is {2}.'.format(x_len, num, pdb_x_len))

        # The space at x side between two molecules
        space = (x_len - (num + 0.5) * pdb_x_len) / (num - 0.5)
        if space < 1.5:
            echo_warning('The space %s between two molecules is less than 1.5. It maybe take a long time to pack.' % space)

        begin_coord = init_coord if even else init_coord + 0.5 * (space + pdb_x_len)
        end_coord = begin_coord + num * (pdb_x_len + space)
        return (begin_coord, end_coord)

    def __mem_loose_pack(self):
        '''Packing molecules random at xy area, but restrict some atoms at z direction.'''
        echo_note('Packing {0[0]}*{0[1]} molecules in the box of {1[0]} {1[1]} {1[2]} {1[3]} {1[4]} {1[5]}.'.format(self.pack['align_num'], self.pack['box']))

        self.results.append('')
        x1, x2 = self.pack['box'][0], self.pack['box'][3]
        y1, y2 = self.pack['box'][1], self.pack['box'][4]
        self.results[0] += 'structure %s\n' % self.pack['pdb']
        # 添加殘基排序选项以使膜堆砌时残基序号递增
        self.results[0] += '\tresnumbers 2 \n'
        self.results[0] += '\tnumber ' + str(self.pack['align_num'][0] * self.pack['align_num'][1]) + '\n'
        self.results[0] += '\tinside box'
        z1 = self.pack['box'][2]
        z2 = self.pack['box'][5]
        for m in x1, y1, z1, x2, y2, z2:
                self.results[0] += ' ' + trans_num(m)
        self.results[0] += '\n'
        self.results[0] += self.__set_panel(x1, y1, z1, x2, y2, z2)
        self.results[0] += 'end structure\n\n'

    def __random_pack(self):
        '''Packing molecules in an area, no atom restriction.'''
        echo_note('Packing {0} molecules in the box of {1[0]} {1[1]} {1[2]} {1[3]} {1[4]} {1[5]}.'.format(self.pack['num'], self.pack['box']))
        self.results.append('structure {0}\n\tnumber {1}\n\tinside box {2[0]} {2[1]} {2[2]} {2[3]} {2[4]} {2[5]}\nend structure\n\n'.format(self.pack['pdb'], self.pack['num'], self.pack['box']))


class ReadConfig(configparser.ConfigParser):
    '''Parse ini configuration file.

        The dict structure is as below:
            {'pack1':
                {'pdb':'m',
                 'box':(0, 0, 0, 10, 10, 10),
                 'fixed': True/False
                 'membrane': True/False,
                 'align_num': (2,3)
                 'takeout': 3
                 'num': 10,
                 'over_x':{10:[23,45],...},
                 'below_x':{10:[23,45],...},
                 'over_x':{10:[23,45],...},
                 'below_x':{10:[23,45],...},
                 'over_x':{10:[23,45],...},
                 'below_x':{10:[23,45],...}
                 },
            ...}
            pdb: pdb file name.
            box: The six size values of box. e.g. [0,0,0,10,10,10]
            fixed: Wheather to fix this pdb，It can be used to fix a pdb as its origin position but not to pack it. When it is true, only the option pdb are available. The number of pdb will be set to 1.
            membrane: Wheather to construct a membrane. The membrane will paralell to the z axel.
            align_num: e.g. (5,4), if membrane option is true, the align_num is the number of x and y axles, and the num option is unavailable. The diagram is as below:
                  y     .   .   .   .
                  y       .   .   .   .
                  y     .   .   .   .
                  y       .   .   .   .
                  y     .   .   .   .
                        x   x   x   x
                  This is a orthorhombic which x*y=4*5, totally 20 molecules.
            takeout: subtract several molecules from above packing. It can pack the molecules with a number is not regular, such as number between 6*5 and 6*4. This option is only for membranePack.
            num: The number of molecules, only when membrane option is false.
            over_x: Make atom above plane paralelling x axle. The value is such as {x coordinate:[atom index, ...], ...}， e.g. make the indices 23,45 atoms above the plane x=10 can be writen as {10:[23,45]}.
            below_x: similar to above option but below plane paralelling x axle.
            over_y: ...
            below_y: ...
            over_z: ...
            below_z: ...
        '''
    def __init__(self, fin):
        super(ReadConfig, self).__init__(inline_comment_prefixes=('#', ';'))
        self.read_file(fin)
        self.packs = OrderedDict()
        self.__parser()

    def __parser(self):
        '''Parser ini file to get packing options'''
        # The number of kinds of molecules
        self.pack_num = self.getint('packmol', 'pack_num')
        # The final file name of pdb
        self.out_pdb = self.get('packmol', 'out_pdb')

        # Read every [pack*] section
        for i in range(self.pack_num):
            sec = 'pack' + str(i)
            secs = {}
            # pdb file name
            secs['pdb'] = self.get(sec, 'pdb')
            # Pack type: membrane, fixed, random
            secs['pack_type'] = self.get(sec, 'pack_type')

            if secs['pack_type'] == 'membrane' or secs['pack_type'] == 'random':
                # box size
                secs['box'] = [float(k) for k in self.get(sec, 'box', fallback='0 0 0 0 0 0').split()]
                # Restriction of atom
                for l in ['over', 'below']:
                    for m in ['x', 'y', 'z']:
                        panel_type = l + '_' + m
                        panel_temp = self.get(sec, panel_type, fallback='')
                        if len(panel_temp) > 0:
                            secs[panel_type] = eval(panel_temp)
                        else:
                            secs[panel_type] = {}
                if secs['pack_type'] == 'membrane':
                    # packing method
                    secs['method'] = self.get(sec, 'method', fallback='copy')
                    # The number in x and y direction.
                    secs['align_num'] = [int(k) for k in self[sec]['align_num'].split()]
                    # Subtract some molecules
                    secs['takeout'] = self.getint(sec, 'takeout', fallback=0)
                else:
                    secs['num'] = self.getint(sec, 'num', fallback=0)
            elif secs['pack_type'] == 'fixed':
                # molecule number
                secs['num'] = 1
                # offset
                secs['offset'] = [float(k) for k in self.get(sec, 'offset', fallback='0. 0. 0. 0. 0. 0.').split()]
            else:
                echo_error('The pack_type must be one of membrane, fixed, random in %s' % sec)
                exit(1)

            self.packs[sec] = secs

    def __valid_check(self):
        '''Check the validity of options.'''
        if len(self.out_pdb) == 0:
            echo_error("You must assign the output pdb file name in [packmol].")
            exit(1)

        # Read every [pack*] section
        for sec, options in self.packs.items():
            # The number of molecules
            if options['pack_type'] == 'membrane':
                if not isinstance(options['align_num'][0] * options['align_num'][1], int):
                    echo_error('The align_num of molecules must be integer in [%s].' % sec)
                    exit(1)
                if not options['method'] in ('copy', 'general', 'line', 'loose'):
                    echo_error('The method of packing membrane must be one of copy, general, line, loose in ' % sec)
                    exit(1)
            elif not isinstance(options['num'], int):
                echo_error('The number of molecules must be integer in [{}].'.format(sec))
                exit(1)

            if options['pack_type'] == 'membrane' or options['pack_type'] == 'random':
                if len(options['box']) != 6:
                    echo_error('The box size should have six values in [%s]' % sec)
                    exit(1)
            elif options['pack_type'] == 'fixed':
                if len(options['offset']) != 6:
                    echo_error('The offset should have six values in [%s]' % sec)
                    exit(1)

            panel_types = ['over_x', 'below_x', 'over_y', 'below_y', 'over_z', 'below_z']
            for option in panel_types:
                if option in self.packs.keys():
                    value = self.packs[option]
                    # Check coordinate and atom indices
                    for coor, atoms in value.items():
                        if not isinstance(coor, (int, float)):
                            echo_error('The coordinate must be integer or float for {} option in {}'.format(option, sec))
                            exit(1)
                        if len(atoms) == 0:
                            echo_error('The number of atoms must be bigger than 0 in {} option in {}'.format(option, sec))
                            exit(1)
                        for atom in atoms:
                            if not isinstance(atom, int) or atom < 0:
                                echo_error('The index of atom must be integer.')
                                exit(1)


def main():
    arg_parser = argparse.ArgumentParser(description='Read ini configure file and packing molecules.')
    arg_parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True, help='Ini configure file.')
    arg_parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True, help='Output inp file.')
    arg_parser.add_argument('-l', '--loop', type=int, action='store', default=100000, help='nloop value, which is max interate number.')
    arg_parser.add_argument('-t', '--tolerance', type=float, action='store', default=2.0, help='tolerance, the mininum distance between two molecules. Decreasing this value can packing quicker but maybe make molecules too tense, default is 2.0. It should be no less than 1.5.')
    arg_parser.add_argument('--parallel', type=int, action='store', help='Whether or not to use paralell version.')
    arg_parser.add_argument('--nopause', action='store_true', help='No pausing confirmation before packing')

    args = arg_parser.parse_args()

    f = args.output

    if args.loop:
        loop = args.loop
    if args.tolerance:
        tolerance = args.tolerance

    config_reader = ReadConfig(args.input)
    f.write('tolerance {}\noutput {}\nfiletype pdb\nnloop {}\n\n'.format(tolerance, config_reader.out_pdb, loop))
    packing_files = []
    for pack in config_reader.packs.values():
        mol_pack = MolPack(pack)
        results = mol_pack.packmol()
        if len(results) == 0:
            echo_error("No results generated. Please check the configure file.")
            exit(1)
        f.write(results.pop(-1))
        for index, result in enumerate(results):
            now = str(time.time()).replace('.', '')
            former_inp = '/tmp/mid' + now + '.inp'
            with open(former_inp, 'w') as ff:
                ff.write('tolerance {}\nfiletype pdb\nnloop {}\n'.format(tolerance, loop))
                ff.write(result)
            packing_files.append(former_inp)
    packing_files.append(f.name)
    f.close()
    for inp in packing_files:
        if args.parallel:
            ppackmol(inp, args.parallel, args.nopause)
        else:
            packmol(inp, nopause=args.nopause)


if __name__ == '__main__':
    main()
