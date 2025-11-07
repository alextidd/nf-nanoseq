import re

class GInterval:
    def __init__(self, chrr, beg, end):
        self.chr = chrr
        if end < beg:
            raise ValueError("Interval %s: %s - %s is invalid!" % (chrr, beg, end))
        self.beg = beg
        self.end = end
        self.l = end - beg + 1

    def convert2DSAInput(self):
        return GInterval(self.chr, self.beg-1, self.end-1)

    def __repr__(self):
        return repr((self.chr, self.beg, self.end))

    def __str__(self):
        return "%s:%s-%s" % (self.chr, self.beg, self.end)

    def __contains__(self, other):
        if self.chr != other.chr:
            return False
        if other.end < self.beg or self.end < other.beg:
            return False
        return True

    def __sub__(self, other):
        if self.chr != other.chr:
            return [self]
        if other.beg > self.beg and other.beg <= self.end:
            if other.end < self.end:
                return [GInterval(self.chr, self.beg, other.beg-1), GInterval(self.chr, other.end+1, self.end)]
            if other.end >= self.end:
                return [GInterval(self.chr, self.beg, other.beg - 1)]
        if other.beg <= self.beg:
            if other.end < self.end and other.end >= self.beg:
                return [GInterval(self.chr, other.end + 1, self.end)]
            if other.end >= self.end:
                return []
        return [self]

    def __add__(self, other):
        if other.chr != self.chr:
            return [self, other]
        if other.end + 1 == self.beg or self.end + 1 == other.beg:
            return [GInterval(self.chr, min([self.beg, other.beg]), max([self.end, other.end]))]
        if other.end < self.beg or self.end < other.beg:
            return [self, other]
        if self.beg <= other.end or other.beg <= self.end:
            return [GInterval(self.chr, min([self.beg, other.beg]), max([self.end, other.end]))]

    @staticmethod
    def chr_isless(x, y):
        if x == y:
            return False
        xc = re.sub('^chr', '', x)
        yc = re.sub('^chr', '', y)
        if xc.isdigit():
            if yc.isdigit():
                return int(xc) < int(yc)
            else:
                return True
        elif yc.isdigit():
            return False
        elif (xc == "M" or xc == "MT") and (yc == "X" or yc == "Y"):
            return False
        elif (xc == "X" or xc == "Y") and (yc == "M" or yc == "MT"):
            return True
        elif xc == "X":
            return True
        elif yc == "X":
            return False
        elif xc == "Y":
            return True
        elif yc == "Y":
            return False
        elif xc == "M" or xc == "MT":
            return True
        elif yc == "M" or yc == "MT":
            return False
        else:
            return xc < yc

    def __eq__(self, other):
        if self.chr == other.chr:
            if self.beg == other.beg:
                if self.end == other.end:
                    return True
        return False

    def __ne__(self, other):
        if self.chr == other.chr:
            if self.beg == other.beg:
                return False
        return True

    def __gt__(self, other):
        if GInterval.chr_isless(other.chr, self.chr):
            return True
        if GInterval.chr_isless(self.chr, other.chr):
            return False
        if self.beg > other.beg:
            return True
        else:
            return False

    def __ge__(self, other):
        if GInterval.chr_isless(other.chr, self.chr):
            return True
        if GInterval.chr_isless(self.chr, other.chr):
            return False
        if self.beg >= other.beg:
            return True
        else:
            return False

    def __lt__(self, other):
        if GInterval.chr_isless(self.chr, other.chr):
            return True
        if GInterval.chr_isless(other.chr, self.chr):
            return False
        if self.beg < other.beg:
            return True
        else:
            return False

    def __le__(self, other):
        if GInterval.chr_isless(self.chr, other.chr):
            return True
        if GInterval.chr_isless(other.chr, self.chr):
            return False
        if self.beg <= other.beg:
            return True
        else:
            return False