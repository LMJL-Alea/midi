# In MIDI, we use the `constants.R` file to store default values for physical
# constants and compartment parameters. This file is sourced by the
# `compartment.R` file to set default values for the compartment parameters.

# Units:
# - um
# - ms
# - mT

gyromagnetic_ratio <- function() {
  267.5987 # rad ms^-1 mT^-1
}

default_free_diffusivity <- function() {
  3.0 # um^2 ms^-1
}

default_sphere_radius <- function() {
  15 # um
}

default_cylinder_radius <- function() {
  0.5 # um
}
