"""
##############################################################################
#                        KITE | Release                                      #
#                                                                            #
#  Developed by: Henrique P. Veiga, Joao V. Lopes,                           #
#  Kevin J. U. Vidarte, Tatiana G. Rappoport, Aires Ferreira, 2025           #
#                                                                            #
##############################################################################
"""

class Vertex:
    def __init__(self, n_, stream_ = []):
        """Tn Operator"""
        self.moment = n_
        self.stream = stream_
    def add_operator(self, coef_, str_):
        """coef_ * Operator -> str_"""
        self.stream.append([coef_, str_])
