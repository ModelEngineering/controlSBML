"""Encoding of options as a dictionary"""


class Options(dict):
    # Class to manage options

    def __init__(self, dct=None):
        super().__init__(dct)

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
        is_default = False
        if key in self.keys():
            is_default = self[key] is None
        else:
            is_default = True
        if is_default:
            self[key] = default

    def parse(self, dcts):
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
        all_opts = set([])
        for dct in dcts:
            all_opts = all_opts.union(dct.keys())
        unknownOptions = set(self.keys()).difference(all_opts)
        if len(unknownOptions) > 0:
            raise ValueError("Unknown options: %s" % str(unknownOptions))
        #
        opts_lst = []
        for dct in dcts:
            opts = Options({k: self[k] if k in self.keys() else v
                  for k, v in dct.items()})
            opts_lst.append(opts)
        return opts_lst
