from collections import namedtuple

# Segment ###########################################
_s_field_names = 'left right node children time population'
_s_values = (0.0, 1.0, None, (), -0.0, 0)
default_segment = dict(zip(_s_field_names.split(' '),
                              _s_values))
Segment = namedtuple('Segment', _s_field_names)
Segment.__new__.__defaults__ = _s_values

def segment_str(self):
    parent = self.node
    if parent is None:
        parent = "ancestor"
    out = "[{sl:.{digits}f}, {sr:.{digits}f}) {sn:{nodedigits}} <- {sc} at {st:.{digits}f} in {sp}".format(
        sl=self.left,
        sr=self.right,
        sn=parent,
        sc=self.children,
        st=self.time,
        sp=self.population,
        digits=2,
        nodedigits=8,
    )
    return out

Segment.__str__ = segment_str
Segment.__repr__ = Segment.__str__


# Chromosome ###########################################
_dc_field_names = 'identity breakpoint lparent rparent left_end right_end'
_dc_values = (None, None, None, 0.0, 1.0)
default_chromosome = dict(zip(_dc_field_names.split(' ')[1:],
                              _dc_values))
Chromosome = namedtuple('Chromosome', _dc_field_names)
Chromosome.__new__.__defaults__ = _dc_values
Chromosome.__doc__ = \
    '''a chromosome is produced by recombination between two parents

        Args:
            identity: integer ID for this chromosome
            breakpoint: separates the right from left parts
            lparent: integer ID for the parent of the left part
            rparent: integer ID for the parent of the left part

        Optional args:
            left_end
            right_end -- NOTE no error checking to make sure bp is between these
    '''


def chromosome_str(self):
    SIZE = 4
    if self.breakpoint is not None \
            and self.rparent is not None \
            and self.lparent is not None:
        bp = "{breakpoint:.{bpdigits}f})(".format(breakpoint=self.breakpoint,
                                                  bpdigits=2)
        info = "{lp:<{ids}}~ {bp} ~{rp:>{ids}}".format(lp=self.lparent,
                                                       rp=self.rparent,
                                                       bp=bp,
                                                       ids=SIZE)
    elif self.breakpoint is None \
            and self.rparent is None \
            and self.lparent is None:
        info = "~ ancestor ~"
    else:
        raise ValueError("lparent, rparent and breakpoint",
                         "must all be None or not None")

    output = "{id:{idsize}}: {le:.{epdigits}f}.[" + info + "]{re:.{epdigits}f}."
    return output.format(
        id=self.identity,
        le=self.left_end,
        re=self.right_end,
        epdigits=0,
        idsize=SIZE
    )

Chromosome.__str__ = chromosome_str
Chromosome.__repr__ = Chromosome.__str__
