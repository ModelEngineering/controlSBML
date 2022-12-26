"""Encoding of options as a dictionary"""

"""
TO DO
1. Construct with default_dcts, list of defaults
2. self.default_dct is the union of default_dcts
3. set checks the default value, and reassigns if at default
"""


class Options(dict):
    # Class to manage options

    def __init__(self, dct, default_dcts):
        if dct is None:
            dct = {}
        super().__init__(dct)
        # Validate the inputs
        if not isinstance(dct, dict):
            raise ValueError("First argument must be a dict.")
        if not isinstance(default_dcts, list):
            raise ValueError("Second argument must be a list of dict.")
        else:
            for default_dct in default_dcts:
                if not isinstance(default_dct, dict):
                    raise ValueError("Second argument must be a list of dict.")
        # Create the set of all defaults
        if default_dcts is None:
            default_dcts = []
        self.default_dcts = default_dcts  # Options  by dict
        self.all_default_dct = {}
        for default_dct in default_dcts:
            self.all_default_dct.update(default_dct)

    def __repr__(self):
        return str({k: v for k, v in self.items()})

    def set(self, key, override=None, default=None):
        """
        Sets the value of an option in an option dictionary.

        Parameters
        ----------
        key: str
        opt_dct: dict
            key: option name
            value: value of the option
        override: obj
            value of setting if not None
        default: obj
            value of setting if keyword is absent or None
        """
        if override is not None:
            self[key] = override
            return
        #
        is_default = False
        if key in self.keys():
            is_default = self[key] == self.all_default_dct[key]
        else:
            is_default = True
        if is_default:
            self[key] = default

    def parse(self):
        """
        Parses options in to lists coresponding to the dictionaries provided.
        All options from each dictionary are included in the results.

        Parameters
        ----------
        kwargs: dict
        dcts: list-dict
            dictionaries to be parsed
        
        Returns
        -------
        list-opts
        """
        # Validate that correct options are present
        unknownOptions = set(self.keys()).difference(self.all_default_dct.keys())
        if len(unknownOptions) > 0:
            raise ValueError("Unknown options: %s" % str(unknownOptions))
        #
        opts_lst = []
        for dct in self.default_dcts:
            opts = Options({k: self[k] if k in self.keys() else v
                  for k, v in dct.items()}, default_dcts=[dct])
            opts_lst.append(opts)
        return opts_lst
