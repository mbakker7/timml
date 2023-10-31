Line-doublets
=============

Line-doublets are used to simulate leaky and impermeable walls. The flux normal to an
impermeable wall is zero. The flux through a leaky wall (the normal component
:math:`q_{n}`) is defined as

    .. math::
        q_n = (h^- - h^+)/c
        
where :math:`h^-` and :math:`h^-` are the heads on the minus and plus sides of the
wall, and :math:`c` is the resistance against flow through the walll. An impermeable
wall is equivalent to a leaky wall with a resistance that is equal to infinity.

1. :class:`~timml.linedoublet.ImpLineDoublet` is used to simulate one straight
   impermeable wall

2. :class:`~timml.linedoublet.ImpLineDoubletString` is an impermeable wall represented
   by a poly line of straight segments

3. :class:`~timml.linedoublet.LeakyLineDoublet` is used to simulate one straight leaky
   wall

4. :class:`~timml.linedoublet.LeakyLineDoubletString` is a leaky wall represented by a
   poly line of straight segments