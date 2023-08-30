


class PDB_Object():
    def __init__(self,file):
        self.file = file
        self.name = file.name
        name, ext = file.name.split('.')
        self.id = name[:4]
        self.no_as_file = dir + '/' + self.id + '.no_info'
        self.no_solv_file = dir + '/' + self.id + '.no_solv'

    def get_residues(self):
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(self.id, self.file)
        self.residues = [r for r in struct.get_residues()]

    def get_rcsb_tree(self):
        ''' Return the xml tree of the rscb page.
        '''
        parser = etree.HTMLParser()
        base_rcsb_url = 'https://www.rcsb.org/structure/'
        url = base_rcsb_url + self.id
        print(url)
        page = requests.get(url)
        html = page.content.decode("utf-8")
        return etree.parse(StringIO(html), parser=parser)

    def check_for_uniprot_id(self):
        ''' Resolve the identity of the uniprot object linked
        to the pdb model by extracting the name from the uniprot
        url that is linked in the rscb page for the model.
        '''
        tree = self.get_rcsb_tree()
        refs = tree.xpath("//a")
        links = [link.get('href', '') for link in refs]
        uniprot_links = [l for l in links if 'uniprot.org' in l]
        if len(uniprot_links) == 0:
            return False
        # Ex: 'https://www.uniprot.org/uniprot/P00942' --> 'P00942'
        self.uniprot_id = uniprot_links[0].split('/')[-1]
        return True

    def fetch_uniprot_text(self):
        ''' Get the raw text of the uniprot text file corresponding
        with the pdb model.
        '''
        uniprot_base_url = 'https://rest.uniprot.org/uniprotkb/'
        uniprot_id = self.uniprot_id
        url = uniprot_base_url + uniprot_id + '.txt'
        page = requests.get(url)
        return page.text

    def resolve_res_ids(self):
        ''' Extract the residue ID numbers of the active site
        residues from the uniprot text file for the pdb model.
        '''
        # Text file from uniprot detailing model features
        raw_text = self.fetch_uniprot_text()
        # Split text into list of lines
        text = raw_text.split('\n')
        # Split each line into list of words
        parsed = (parse.text_to_list(line) for line in text)
        # Isolate lines containing active site info
        act_site_lines = (line for line in parsed if 'ACT_SITE' in line)
        # Residue IDs are the last words in the lines
        return tuple(int(line[-1]) for line in act_site_lines)

    def get_active_site_residues(self):
        ''' Get residues belonging to the active site.
        '''
        act_site_ids = self.resolve_res_ids()
        self.as_res = [r for r in self.residues if r.id[1] in act_site_ids]
        self.as_res = [r for r in self.as_res if not r.resname in ('HOH','SOL')]

    def has_active_site_information(self):
        return len(self.as_res) > 0

    def get_outfile_path(self):
        ''' Return path of output file to make.
        '''
        outfile_template = '{}_wireInfo_cutoff{}_new.csv'
        outfile = outfile_template.format(self.id, int(CUTOFF))
        return os.path.join(dir, outfile)

    def make_no_as_file(self):
        ''' Make a file to indicate that the model has no
        active site information'''
        os.system('touch {}'.format(self.no_as_file))

    def make_no_solv_file(self):
        ''' Make a file to indicate that the model has no
        active site information'''
        os.system('touch {}'.format(self.no_solv_file))

    def has_as_file(self):
        ''' Check if file exists to indicate that no AS info can
        be found.
        '''
        return os.path.exists(self.no_as_file )


#### After running once, taking outputs to check for direct AS-bulk conn ####

mega = {'Enzyme Atoms':[],
        'AS Atoms': [],
        'Wire Length': [],
        'Water Near enz_coords': [],
        'N Solvent': [],
        'N Solvent Contacted': [],
        'Class': [],
        'Filename': []}


def pdb_has_no_contacted_solvent(pdb):



    pdb.outfile = pdb.get_outfile_path()
    if os.path.exists():
        with open(pdb.outifle, 'r') as f:
            lines = [l for l in f]
        del lines[0]
        lines = [l.replace('\n','') for l in lines]
        lines = [l.split(',') for l in lines]
