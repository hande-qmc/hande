Old (removed) functionality
===========================

Unused and **not useful** functionality is occasionally removed from HANDE, in
order to remove the maintenance burden for code that really has no benefit.  In
general, keeping failed experiments in the codebase is not helpful to
developers (more work) and users (not obvious if an option should or should not
be used).  When it transpires that something falls into the category, we may
hence remove it and detail it below.  If you are interested in resurrecting
this functionality, please dig through the git history and/or speak to
a developer.

folded-spectrum FCIQMC
    The folded-spectrum approach allows, in principle, access to excited states
    in FCIQMC via using the Hamiltonian :math:`(H-\epsilon)^2`, where
    :math:`\epsilon` is an energy offset.  It emerged in practice to be very
    painful/impossible to converge to excited states for systems beyond the
    reach of conventional FCI.
defining an initiator determinant via a complete active space
    Originally the initiator space was defined by a population threshold and
    a complete active space (CAS).  It turns out that it is simpler to allow
    the initiator space to emerge naturally just through the population
    threshold (as used in later studies), whereas defining a CAS that is small
    but effective is not easy in large systems.  Furthermore, using just
    a population threshold makes the initiator approximation easier to extend
    to other algorithms (i.e. CCMC and DMQMC).
