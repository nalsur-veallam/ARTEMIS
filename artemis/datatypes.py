import json
import numpy as np


class MAP:
    """
    Class to store MIE matrix information.
    """

    def __init__(self, filename=None, data=None):
        """
        Initialize Map instance. Optionally a file path can be provided to directly read the file.

        Args:
            filename (str): Path to the file containing MIE matrix
        """

        self.errors = []
        self.map_ = None
        self.NResidues = None
        self.names = None
        self.real_numbers = None

        if filename is not None:

            try:
                with open(filename) as json_file:
                    data = json.load(json_file)
            except:
                self.errors.append(ValueError(f"Can't read data from file {filename}."))

            if self.is_ok():
                self.read_data(data, filename)
                if self.exist():
                    self.check()

        elif data is not None:
            self.read_data(data)
            if self.exist():
                self.check()


    def read_data(self, data, filename=None) -> None:

        try:
            self.map_ = np.array(data['map'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'map' key."))

        try:
            self.NResidues = int(data['NResidues'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'NResidues' key."))

        try:
            self.names = np.array(data['names'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'names' key."))

        try:
            self.real_numbers = np.array(data['real_numbers'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'real_numbers' key."))

    def check(self) -> None:

        if len(self.names) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the names array does not match the number of residues."))

        if self.map_.shape[0] != self.NResidues or self.map_.shape[1] != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the MIE matrix does not match the number of residues."))

        if len(self.real_numbers) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the real_numbers array does not match the number of residues."))

    def is_ok(self) -> bool:
        if self.errors == []:
            return True
        else:
            return False

    def interrupt(self):
        if not self.is_ok():
            for error in self.errors:
                raise error

    def exist(self) -> bool:

        if self.names is None:
            return False
        elif self.map_ is None:
            return False
        elif self.NResidues is None:
            return False
        elif self.real_numbers is None:
            return False
        else:
            return True

    def write(self, out_path) -> None:
        if self.exist() and self.is_ok():
            data = {}
            data['ftype'] = 'map'
            data['names'] = self.names.tolist()
            data['NResidues'] = self.NResidues
            data['map'] = self.map_.tolist()
            data['real_numbers'] = self.real_numbers.tolist()
            try:
                with open(out_path, 'w') as outfile:
                    json.dump(data, outfile)
                print("File",out_path + " created\n")
            except:
                raise ValueError(f"Error writing file {out_path}.")


class GROUPS:
    """
    Class to store groups information.
    """

    def __init__(self, filename=None, data=None):
        """
        Initialize GROUPS instance. Optionally a file path can be provided to directly read the file.

        Args:
            filename (str): Path to the file containing MIE matrix
        """

        self.errors = []

        if filename is not None:

            try:
                with open(filename) as json_file:
                    data = json.load(json_file)
            except:
                self.errors.append(ValueError(f"Can't read data from file", filename, "."))

        self.data = data

    def exist(self) -> bool:

        if self.data is None:
            return False
        else:
            return True

    def is_ok(self) -> bool:
        if self.errors == []:
            return True
        else:
            return False

    def interrupt(self):
        if not self.is_ok():
            for error in self.errors:
                raise error

    def write(self, out_path) -> None:
        if self.exist() and self.is_ok():
            try:
                with open(out_path, 'w') as outfile:
                    json.dump(self.data, outfile)
                print("File",out_path + " created\n")
            except:
                raise ValueError(f"Error writing file {out_path}.")

class ALLOSTERY:
    """
    Class to store allostery information.
    """

    def __init__(self, filename=None, data=None):
        """
        Initialize Allostery instance. Optionally a file path can be provided to directly read the file.

        Args:
            filename (str): Path to the file containing MIE matrix
        """

        self.errors = []
        self.map_ = None
        self.NResidues = None
        self.names = None
        self.real_numbers = None
        self.active_site = None
        self.allosteric_site = None
        self.intensity = None

        if filename is not None:

            try:
                with open(filename) as json_file:
                    data = json.load(json_file)
            except:
                self.errors.append(ValueError(f"Can't read data from file {filename}."))

            if self.is_ok():
                self.read_data(data, filename)
                if self.exist():
                    self.check()

        elif data is not None:
            self.read_data(data)
            if self.exist():
                self.check()


    def read_data(self, data, filename=None) -> None:

        try:
            self.map_ = np.array(data['map'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'map' key."))

        try:
            self.NResidues = int(data['NResidues'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'NResidues' key."))

        try:
            self.names = np.array(data['names'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'names' key."))

        try:
            self.real_numbers = np.array(data['real_numbers'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'real_numbers' key."))

        try:
            self.active_site = np.array(data['active_site'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'active_site' key."))

        try:
            self.allosteric_site = np.array(data['allosteric_site'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'allosteric_site' key. Continues without this data.")
            self.allosteric_site = None

        try:
            self.intensity = np.array(data['intensity'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'intensity' key. Continues without this data.")
            self.intensity = None

    def check(self) -> None:

        if len(self.names) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the names array does not match the number of residues."))

        if self.map_.shape[0] != self.NResidues or self.map_.shape[1] != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the MIE matrix does not match the number of residues."))

        if len(self.real_numbers) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the real_numbers array does not match the number of residues."))

        if max(self.active_site) > self.NResidues or min(self.active_site) < 1 or len(self.active_site) == 0:
            self.errors.append(ValueError(f"The received data is broken: incorrect active site format."))

        if self.allosteric_site is not None and (max(self.allosteric_site) > self.NResidues or min(self.allosteric_site) < 1 or len(self.allosteric_site) == 0):
            self.errors.append(ValueError(f"The received data is broken: incorrect allosteric site format."))

        if self.intensity is not None and min(self.intensity) < 0:
            self.errors.append(ValueError(f"The received data is broken: incorrect intensity data."))

    def is_ok(self) -> bool:
        if self.errors == []:
            return True
        else:
            return False

    def interrupt(self):
        if not self.is_ok():
            for error in self.errors:
                raise error

    def exist(self) -> bool:

        if self.names is None:
            return False
        elif self.map_ is None:
            return False
        elif self.NResidues is None:
            return False
        elif self.real_numbers is None:
            return False
        elif self.active_site is None:
            return False
        else:
            return True

    def from_map_and_groups(self, Map, Groups, active_site_key="active_site", allosteric_site_key="allosteric_site") -> None:
        if not Map.exist() or not Map.is_ok():
            self.errors.append(ValueError(f"The received Map object is broken."))
        if not Groups.exist() or not Groups.is_ok():
            self.errors.append(ValueError(f"The received Groups object is broken."))

        if active_site_key in Groups.data.keys():
            self.active_site = np.array(Groups.data[active_site_key])
            if allosteric_site_key in Groups.data.keys():
                self.allosteric_site = np.array(Groups.data[allosteric_site_key])

            self.map_ = Map.map_
            self.NResidues = Map.NResidues
            self.names = Map.names
            self.real_numbers = Map.real_numbers
            if self.exist():
                self.check()

        else:
            self.errors.append(ValueError(f"The received data is broken: can't find active site group by {active_site_key} key."))


    def write(self, out_path) -> None:
        if self.exist() and self.is_ok():
            data = {}
            data['ftype'] = 'allostery'
            data['names'] = self.names.tolist()
            data['NResidues'] = self.NResidues
            data['map'] = self.map_.tolist()
            data['real_numbers'] = self.real_numbers.tolist()
            data['active_site'] = self.active_site.tolist()
            if self.allosteric_site is not None:
                data['allosteric_site'] = self.allosteric_site.tolist()
            if self.intensity is not None:
                data['intensity'] = self.intensity.tolist()
            try:
                with open(out_path, 'w') as outfile:
                    json.dump(data, outfile)
                print("File",out_path + " created\n")
            except:
                raise ValueError(f"Error writing file {out_path}.")

class CLUSTERS:
    """
    Class to store allostery information.
    """

    def __init__(self, filename=None, data=None, NClusters=None):
        """
        Initialize Allostery instance. Optionally a file path can be provided to directly read the file.

        Args:
            filename (str): Path to the file containing MIE matrix
        """

        self.errors = []
        self.map_ = None
        self.NResidues = None
        self.names = None
        self.real_numbers = None
        self.reference_group = None
        self.restriction_group = None
        self.active_site = None
        self.allosteric_site = None
        self.clustering_labels = None
        self.NClusters = NClusters

        if filename is not None:

            try:
                with open(filename) as json_file:
                    data = json.load(json_file)
            except:
                self.errors.append(ValueError(f"Can't read data from file {filename}."))

            if self.is_ok():
                self.read_data(data, filename)
                if self.exist():
                    self.check()

        elif data is not None:
            self.read_data(data)
            if self.exist():
                self.check()


    def read_data(self, data, filename=None) -> None:

        try:
            self.map_ = np.array(data['map'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'map' key."))

        try:
            self.NResidues = int(data['NResidues'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'NResidues' key."))

        try:
            self.names = np.array(data['names'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'names' key."))

        try:
            self.real_numbers = np.array(data['real_numbers'])
        except:
            self.errors.append(ValueError(f"Can't get data from file {filename} by 'real_numbers' key."))

        try:
            self.active_site = np.array(data['active_site'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'active_site' key. Continues without this data.")
            self.active_site = None

        try:
            self.allosteric_site = np.array(data['allosteric_site'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'allosteric_site' key. Continues without this data.")
            self.allosteric_site = None

        try:
            self.reference_group = np.array(data['reference_group'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'reference_group' key. Continues without this data.")
            self.reference_group = None

        try:
            self.restriction_group = np.array(data['restriction_group'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'restriction_group' key. Continues without this data.")
            self.restriction_group = None

        try:
            self.clustering_labels = np.array(data['clustering_labels'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'clustering_labels' key. Continues without this data.")
            self.clustering_labels = None

        try:
            self.NClusters = int(data['NClusters'])
        except:
            print(f"Caution: Can't get data from file {filename} by 'NClusters' key. Continues without this data.")

    def check(self) -> None:

        if len(self.names) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the names array does not match the number of residues."))

        if self.map_.shape[0] != self.NResidues or self.map_.shape[1] != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the MIE matrix does not match the number of residues."))

        if len(self.real_numbers) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the real_numbers array does not match the number of residues."))

        if  self.restriction_group is not None and (max(self.restriction_group) > self.NResidues or min(self.restriction_group) < 1 or len(self.restriction_group) == 0):
            self.errors.append(ValueError(f"The received data is broken: incorrect restriction group format."))

        if self.reference_group is not None and (max(self.reference_group) > self.NResidues or min(self.reference_group) < 1 or len(self.reference_group) == 0):
            self.errors.append(ValueError(f"The received data is broken: incorrect reference_group format."))

        if self.NClusters is not None and self.NClusters <= 0:
            self.errors.append(ValueError(f"The received data is broken: incorrect NClusters format."))

        if self.clustering_labels is not None and len(self.clustering_labels) != self.NResidues:
            self.errors.append(ValueError(f"The received data is broken: the size of the clustering_labels array does not match the number of residues."))

        if self.active_site is not None and (max(self.active_site) > self.NResidues or min(self.active_site) < 1 or len(self.active_site) == 0):
            self.errors.append(ValueError(f"The received data is broken: incorrect active site format."))

        if self.allosteric_site is not None and (max(self.allosteric_site) > self.NResidues or min(self.allosteric_site) < 1 or len(self.allosteric_site) == 0):
            self.errors.append(ValueError(f"The received data is broken: incorrect allosteric site format."))

    def is_ok(self) -> bool:
        if self.errors == []:
            return True
        else:
            return False

    def interrupt(self):
        if not self.is_ok():
            for error in self.errors:
                raise error

    def exist(self) -> bool:

        if self.names is None:
            return False
        elif self.map_ is None:
            return False
        elif self.NResidues is None:
            return False
        elif self.real_numbers is None:
            return False
        elif self.NClusters is None and self.restriction_group is None:
            return False
        else:
            return True

    def from_map_and_groups(self, Map, Groups, reference_group_key="reference_group", restriction_group_key="restriction_group", active_site_key="active_site", allosteric_site_key="allosteric_site") -> None:
        if not Map.exist() or not Map.is_ok():
            self.errors.append(ValueError(f"The received Map object is broken."))
        if not Groups.exist() or not Groups.is_ok():
            self.errors.append(ValueError(f"The received Groups object is broken."))

        if reference_group_key in Groups.data.keys():
            self.reference_group = np.array(Groups.data[reference_group_key])
        if restriction_group_key in Groups.data.keys():
            self.restriction_group = np.array(Groups.data[restriction_group_key])

        if active_site_key in Groups.data.keys():
            self.active_site = np.array(Groups.data[active_site_key])
        if allosteric_site_key in Groups.data.keys():
            self.allosteric_site = np.array(Groups.data[allosteric_site_key])

        self.map_ = Map.map_
        self.NResidues = Map.NResidues
        self.names = Map.names
        self.real_numbers = Map.real_numbers
        if self.exist():
            self.check()

    def read_groups(self, Groups, reference_group_key="reference_group", restriction_group_key="restriction_group", active_site_key="active_site", allosteric_site_key="allosteric_site") -> None:
        if reference_group_key in Groups.data.keys():
            self.reference_group = np.array(Groups.data[reference_group_key])
        if restriction_group_key in Groups.data.keys():
            self.restriction_group = np.array(Groups.data[restriction_group_key])

        if active_site_key in Groups.data.keys():
            self.active_site = np.array(Groups.data[active_site_key])
        if allosteric_site_key in Groups.data.keys():
            self.allosteric_site = np.array(Groups.data[allosteric_site_key])

        if self.exist():
            self.check()


    def write(self, out_path) -> None:
        if self.exist() and self.is_ok():
            data = {}
            data['ftype'] = 'clusters'
            data['names'] = self.names.tolist()
            data['NResidues'] = self.NResidues
            if self.NClusters is not None:
                data['NClusters'] = self.NClusters
            data['map'] = self.map_.tolist()
            data['real_numbers'] = self.real_numbers.tolist()
            if self.active_site is not None:
                data['active_site'] = self.active_site.tolist()
            if self.allosteric_site is not None:
                data['allosteric_site'] = self.allosteric_site.tolist()
            if self.reference_group is not None:
                data['reference_group'] = self.reference_group.tolist()
            if self.restriction_group is not None:
                data['restriction_group'] = self.restriction_group.tolist()
            if self.clustering_labels is not None:
                data['clustering_labels'] = self.clustering_labels.tolist()
            try:
                with open(out_path, 'w') as outfile:
                    json.dump(data, outfile)
                print("File",out_path + " created\n")
            except:
                raise ValueError(f"Error writing file {out_path}.")
