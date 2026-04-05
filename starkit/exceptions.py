class StarKITException(Exception):
    pass


class InvalidInputFileFormat(StarKITException):
    pass


class MissingGFFError(StarKITException):
    pass


class NoProteinsError(StarKITException):
    pass
