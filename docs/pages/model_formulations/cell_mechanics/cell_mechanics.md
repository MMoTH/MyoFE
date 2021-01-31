---
page-title: Cell Level Modeling
title: Cell Mechanics
parent: Model Formulations
has_children: true
nav_order: 1
---

<div class="notice--info">
  <h4>Message</h4>
  <p>This page is under  construction.</p>
</div>

Each gauss point within a finite element is treated as a half-sarcomere. The half-sarcomere is modeled via an embedded Python version of [MyoSim](http://www.myosim.org). MyoSim models half-sarcomeres with cross-bridge distribution techniques. In this implementation, MyoSim is uses information about half-sarcomere length and intracellular calcium concentration [Ca<sup>2+</sup>] to calculate the active stress generated by a representative half-sarcomere. This stress is then scaled according to the density of half-sarcomeres in tissue and assigned to the local fiber direction.  

### Cross-brige Kinetics
This will go more in depth in the equations/cross-bridge distribution techniques.  
[Contraction Models](../cell_mechanics/contraction_models/contraction_models.md)

### Cell Ion Modeling
For now, we are just calculating the calcium concentration. Expand on this.