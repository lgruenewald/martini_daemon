from .__token import Token
import re

class TokenParseError(Exception):
    def __init__(self, token: Token, message: str):
        self.token = token
        self.message = message


class TokenList:
    int_pat = re.compile("[-+]?[0-9]+")
    float_pat = re.compile("[-+]?[0-9]+(\\.[0-9]*)?([eE][-+]?[0-9]+)?")
    word_pat = re.compile("[a-zA-Z0-9_.]+")
    pattern_pat = re.compile("[a-zA-Z0-9_?!*{}]+")
    pair_pat = re.compile(r"[0-9]+:[a-zA-Z0-9_.]+")

    def __init__(self, line: str):
        self.__line = line
        self.__tokens = tokens

    def get_line(self) -> str:
        return self.__line

    def unwrap(self, index, type_filter, default="default placeholder"):
        """
        Given a TokenList try to index it and convert to a usable value
        based on type_filter. If the index is out of range, a default
        value can be specified in place.

        Possible filters:
        int - returns int type, no processing
        float - returns float type, no processing
        positive - float, but raises an exception if 0 or smaller
        index - int, but subtracts 1, exception if 0 or smaller
        degree - degree to radian and makes sure it's in the range 0 to 2pi
        word - a string with only alphanumerics and no whitespace in it
        pattern - fnmatch pattern, converts {} to [] so {} can be used for sets
        pair - a tuple of an index (which reactant) and a word (which graph atom)
        """
        if len(self.__tokens) <= index:
            if default == "default placeholder":
                # hack so "None" can also be used as a default value
                raise TokenParseError(
                    self.__tokens[-1],
                    f"Not enough tokens, expected token at index {index}."
                )
            else:
                return default

        tok = self.__tokens[index].content
        match type_filter:
            case "int":
                if int_pat.match(tok):
                    return int(tok)
            case "float":
                if float_pat.match(tok):
                    return float(tok)
            case "positive":
                if float_pat.match(tok):
                    if float(tok) <= 0.:
                        raise TokenParseError(
                            tokens[index],
                            "Expected a positive non-zero real number."
                        )
                    return float(tok)
            case "index":
                if int_pat.match(tok):
                    if int(tok) <= 0:
                        raise TokenParseError(
                            tokens[index],
                            "Expected index, got an integer 0 or smaller."
                            "Note: indexing in .itp/.top files is usually 1 based."
                        )
                    return int(tok) - 1
            case "degree":
                if float_pat.match(tok):
                    angle = float(tok) * math.pi / 180.0
                    while angle < 0.:
                        angle += math.tau
                    while angle > math.tau:
                        angle -= math.tau
                    return angle
            case "word":
                if word_pat.match(tok):
                    return tok
            case "pattern":
                if pattern_pat.match(tok):
                    return tok.replace("{", "[").replace("}", "]")
            case "pair":
                # specialized index:word construct for [reaction] stuff
                if pair_pat.match(tok):
                    items = tok.split(":")
                    return (int(items[0]) - 1, items[1])

        raise TokenParseError(
            tokens[index],
            f"Expected token of type {type_filter}."
        )
