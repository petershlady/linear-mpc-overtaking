This repository extends a MATLAB implementation of [Safe driving envelopes for path tracking in autonomous vehicles](https://www.sciencedirect.com/science/article/pii/S0967066116300831).

That repository is located [here](https://github.com/petershlady/AutonomousDrivingEnvelopes).

This is used for a course project for the Stanford course AA203 - Optimal Control. We attempted to control an autonomous vehicle overtaking a lead vehicle, possibly in the presence of an oncoming vehicle. Deviation from the path was lightly penalized and soft constraints on stability and obstacle avoidance with much higher costs were included. In the abscense of an oncoming vehicle, the linear MPC formulation (chosen for it's convexity), does well at the overtaking. However, with an oncoming vehicle, it can fail. To preserve linear dynamics, we tried several schemes to account for the fact that the longitudinal speed was variable in the dynamics. However, they are, to a large extent, brittle and less effective. Future work will move to non-linear models.

## Prerequisites

```
MATLAB (Tested with MATLAB R2018b)
CVX (Tested with version 3.0beta)
```
When this was tested with [CVX](http://cvxr.com/cvx/) for MATLAB on Windows 10, it was necessary to use the SeDuMi solver. On MacOS and Ubuntu 16.04, the default SDPT3 solver worker fine. Try them both.

## Authors
This work was performed by Peter Schleede (host of this repo), [Andrew Shoats](https://github.com/ashoats), and [Elliot Weiss](https://github.com/elliotdw).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details