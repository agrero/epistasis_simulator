from epistasiscomponents import constants


class Primogenitor:
    """
The parent gene that will be undergoing mutations in order to create its progeny of Bairn
    """
    def __init__(self, sequence:str, name:str):

        self.sequence = sequence
        self.name = name


    def __str__(self) -> str:
        #can make this a bit more informative, use a few more lines
        #potentially go down the list of variables inputting them as necessary, example below
        #name: {self.name}
        #sequence: {self.sequence}
        return "Parent {} amino acid sequence of {} acids".format(self.name, len(self.sequence))