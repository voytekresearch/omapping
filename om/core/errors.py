"""Custom errors for OM."""

########################################################################################
################################## OM - BASE - ERRORS ##################################
########################################################################################

class UnknownDataSourceError(Exception):
    """An Error indicating data source specification is not understood."""
    pass

class UnknownDataTypeError(Exception):
    """An Error indicating data type specification is not understood."""
    pass

class InconsistentDataError(Exception):
    """An Error indicating there is a fatal inconsistency in data."""
    pass

class DataNotComputedError(Exception):
    """An Error indicating some required data has not been computed."""
    pass
