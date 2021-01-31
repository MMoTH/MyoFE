---
title: Cylinder Twitch
parent: Instruction Files
nav_order: 1
---

<div class="notice--info">
  <h4>Message</h4>
  <p>This page is under  construction.</p>
</div>
```json
{
  "simulation_parameters": {
    "simulation_geometry": ["cylinder"],
    "geometry_options": {
      "base_x": [0.0, "x-coordinate of cylinder base face"],
      "base_y": [0.0, "y-coordinate of cylinder base face"],
      "base_z": [0.0, "z-coordinate of cylinder base face"],
      "base_radius": [1.0],
       "end_x": [10.0, "x-coordinate of cylinder end face"],
       "end_y": [0.0, "y-coordinate of cylinder end face"],
       "end_z":[0.0, "z-coordinate of cylinder end face"],
       "end_radius": [1.0],
       "segments": [20, "number of segments to approximate cylinder"],
       "refinement": [30, "mesh refinement. This combination of segments and refinement yield nodes on y and z axes"]
    },
    "protocol": {
      "simulation_type": ["twitch"],
      "simulation_duration": [150, "ms"],
      "simulation_timestep": [0.5, "ms"]
      },
    "save_cell_output": [1],
    "save_visual_output": [1],
    "output_path": ["./output"]
  },

  "forms_parameters": {

    "passive_law_parameters": {
      "passive_law": ["semi_structural"],
      "c": [4000, "(Pa), isotropic scaling factor"],
      "c2": [500, "(Pa), myofiber stiffness"],
      "c3": [20.0, "(unitless), myofiber exponential factor"],
      "bf": [10.48, "(unitless),Guccione fiber factor"],
      "bt": [3.58, "(unitless), Guccione transverse factor"],
      "bfs": [1.627, "(unitless), Guccione shear factor"],
      "phi_m": [1.0, "(unitless), scaling factor for myofiber passive stress"],
      "phi_g": [1.0, "(unitless), scaling factor for guccione passive stress"]
    }
  },

  "myosim_parameters": {
    "max_rate": [5000,"s^-1"],
    "temperature": [310, "Kelvin"],
    "cb_number_density": [6.96e16, "number of cb's/m^2"],
    "initial_hs_length": [950, "nm"],

    "myofilament_parameters": {
      "kinetic_scheme": ["3state_with_SRX"],
      "num_states": [3],
      "num_attached_states": [1],
      "num_transitions": [4],
      "cb_extensions": [[0.0, 0.0, 4.75642], "power-stroke distance in nm"],
      "state_attached": [[0, 0, 1]],
      "k_cb_multiplier": [[1.0, 1.0, 1.0]],
      "k_cb_pos": [0.001, "N*m^-1"],
      "k_cb_neg": [0.001, "N*m^-1"],
      "alpha":[1.0],
      "k_1": [1.92, "s^-1"],
      "k_force": [4.96e-4, "(N^-1)(m^2)"],
      "k_2": [200.0, "s^-1"],
      "k_3": [272.0, "(nm^-1)(s^-1)"],
      "k_4_0": [508.864648098121948, "s^-1"],
      "k_4_1": [0.189, "nm^-4"],
      "k_cb": [0.001, "N*m^-1"],
      "x_ps": [4.76, "nm"],
      "k_on": [6.53e8, "(M^-1)(s^-1)"],
      "k_off": [200, "s^-1"],
      "k_coop": [6.38],
      "bin_min": [-10, "nm"],
      "bin_max": [10, "nm"],
      "bin_width": [1.0, "nm"],
      "filament_compliance_factor": [0.5],
      "thick_filament_length": [815, "nm"],
      "thin_filament_length": [1120, "nm"],
      "bare_zone_length": [80, "nm"],
      "k_falloff": [0.0024],
      "passive_mode": ["exponential"],
      "passive_exp_sigma": [500],
      "passive_exp_L": [80],
      "passive_l_slack": [950, "nm"],
      "xfiber_fraction": [0, "unitless"]
  }
},

"electrophys_parameters": {

  "cell_ion_parameters": {
    "model": ["two_compartment"],
    "model_inputs": {
      "Ca_content": [1e-3,"Molar"],
      "k_leak": [0.008000800080008],
      "k_act": [1.4472186384349266],
      "k_serca": [80.0]
    }
  },
}
}
```

<a href="../running_a_simulation/running_demo.html" class="btn btn--primary"><< Running a Simulation</a>
<a href="../building_a_mesh/mesh_generation_readme.html" class="btn btn--primary">Building a Mesh >></a>