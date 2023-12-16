import pandas as pd
import controlSBML.constants as cn

COL_CLOSED_LOOP_SYSTEM = "closed_loop_system"
COL_CLOSED_LOOP_SYSTEM_TS = "closed_loop_system_ts"

class History(object):
    # Maintains history of changes to design choices
    def __init__(self, designer, is_history=True):
        self.designer = designer
        self.is_history = is_history
        self._dct = None
        self.clear()

    def __len__(self):
        first = cn.CONTROL_PARAMETERS[0]
        return len(self._dct[first])
    
    def clear(self):
        self._dct = {}
        for name in cn.CONTROL_PARAMETERS:
            self._dct[name] = []
        self._dct[COL_CLOSED_LOOP_SYSTEM] = []
        self._dct[cn.SETPOINT] = []
        self._dct[cn.MSE] = []

    def add(self):
        if not self.is_history:
            return
        for name in cn.CONTROL_PARAMETERS:
            self._dct[name].append(self.designer.__getattribute__(name))
        self._dct[COL_CLOSED_LOOP_SYSTEM].append(self.designer.closed_loop_system)
        self._dct[cn.SETPOINT].append(self.designer.setpoint)
        self._dct[cn.MSE].append(self.designer.residual_mse)

    def undo(self):
        _ = self._dct.pop()

    def report(self):
        """
        Creates a dataframe of the history

        Returns:
            pd.DataFrame
        """
        df = pd.DataFrame(self._dct)
        return df

    def get(self, idx):
        """
        Returns the SISOClosedLoopDesigner at the specified index.

        Args:
            idx: int
        Returns:
            SISOClosedLoopDesigner
        """
        if idx > len(self) - 1:
            raise ValueError("idx must be less than %d" % len(self))
        # Construct entries for the desired history element
        dct = {}
        for name in self._dct.keys():
            dct[name] = self._dct[name][idx]
        designer = SISOClosedLoopDesigner(self.designer.system, self.designer.open_loop_transfer_function,
                                          times=self.designer.times,
                                          setpoint=cn.SETPOINT, is_steady_state=self.designer.is_steady_state,
                                          is_history=self.designer.history.is_history, sign=self.designer.sign)
        for name in cn.CONTROL_PARAMETERS:
            designer.__setattr__(name, dct[name])
        designer.closed_loop_system = dct[COL_CLOSED_LOOP_SYSTEM]
        designer.residual_mse = dct[cn.MSE]
        designer.history.add()
        return designer