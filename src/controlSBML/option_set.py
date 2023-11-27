"""An option set is a set of properties and values. Inherent from this class to incorporate these managed properties.

    Usage:
       class myClass(OptionSet):
            def __init__(self, **kwargs):
                super().__init__(**kwargs) . # Set default values of options

        def method(self, **kwargs):
            option_set = self.getOptionSet(**kwargs)
            actionMethod(option_set.opt1)
 """

class OptionSet(dict):
    # Attributs are also keys

    def __init__(self, **kwargs):
        self.names = list(kwargs.keys())
        for key, value in kwargs.items():
            setattr(self, key, value)
            self[key] = value

    def setOptionSet(self, **kwargs):
        for key, value in kwargs.items():
            self[key] = value
            setattr(self, key, value)

    def _copyOptionSet(self):
        #return OptionSet(**{k: getattr(self, k) for k in self.names})
        option_set = OptionSet(**self)
        return option_set
    
    def __eq__(self, other):
        if not isinstance(other, OptionSet):
            return False
        return self.__dict__ == other.__dict__

    def getOptionSet(self, **kwargs):
        """
        Return a new Option object that overrides as specified by kwargs.

        Args:
            **kwargs: keyword arguments that override the current options
        Returns:
            OptionSet
        """
        new_option_set = self._copyOptionSet()
        new_option_set.setOptionSet(**kwargs)
        return new_option_set