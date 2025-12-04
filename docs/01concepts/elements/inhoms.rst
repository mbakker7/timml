Inhomogeneities
===============

The following polygonal inhomogeneities are available in TimML. Note that there are 
inhomogeneities for both multi-aquifer style models (:class:`~timml.model.ModelMaq`)
and single aquifer models (with multiple layers, :class:`~timml.model.Model3D`).

1. :class:`~timml.inhomogeneity.PolygonInhomMaq`, which is an inhomogeneity consisting
   of a regular sequence of aquifer - leaky layer - aquifer - leaky layer, aquifer, etc.
   The top of the system can be either an aquifer or a leaky layer.

2. :class:`~timml.inhomogeneity.PolygonInhom3D`, which is an inhomogeneity consisting
   of a stack of aquifer layers. The resistance between the aquifer layers is computed as
   the resistance from the middle of one layer to the middle of the next layer. Vertical
   anisotropy can be specified. The system may be bounded on top by a leaky layer.

3. :class:`~timml.inhomogeneity.BuildingPitMaq`, which is an inhomogeneity similar to
   `PolygonInhomMaq` in which impermeable walls can be placed along the edges of the
   inhomogeneity in specified layers.

4. :class:`~timml.inhomogeneity.BuildingPit3D`, which is an inhomogeneity similar to
   `PolygonInhom3D` in which impermeable walls can be placed along the edges of the
   inhomogeneity in specified layers.

5. :class:`~timml.inhomogeneity.LeakyBuildingPitMaq`, which is an inhomogeneity like
   `BuildingPitMaq` but with leaky walls instead of impermeable walls. The resistance of the
   leaky walls can be specified.

6. :class:`~timml.inhomogeneity.LeakyBuildingPit3D`, which is an inhomogeneity like
   `BuildingPit3D` but with leaky walls instead of impermeable walls. The resistance of the
   leaky walls can be specified.