# Wave Breaker

The goal is to simulate proper quantum interactions in a 2D field that can reproduce at least the double slit experiment.


## Simulating Waves
### Problem Statement

Given (electric) wave field 

$$
\psi(\mathbf{r},t): \R^2 \times \mathbb{R} \rightarrow \mathbb{C}
$$

subject to

$$
\nabla^2 \psi = \frac{1}{c^2} \frac{ \partial^2 \psi}{\partial t^2} \ \Leftrightarrow \ c^2 \cdot \nabla^2 \psi = \frac{ \partial^2 \psi}{\partial t^2}
$$

with known initial condidtion $\psi(\mathbf{x},0)$ and $\dot{\psi}(\mathbf{x},0)$ and wave propagation speed $c$, compute in discrete timesteps the evolution of the field

$$
\psi(\mathbf{r},t + \Delta t)
$$

---

Note that the analytical solution to the problem is:

$$
\psi(\mathbf{r},t) = A \cdot e ^{i (\mathbf{k \cdot r} - \omega t)}, \omega = |\mathbf{k}| v
$$

### Analytical Approach

First we have to solve for the second-order spatial partial derivative (Laplacian):

$$
\nabla^2 \psi 
$$

which is equal to the second-order temporal deriviative:

$$
\frac{1}{c^2} \frac{ \partial^2 \psi}{\partial t^2}
$$

by performing a temporal integral of the Laplacian, we get

$$
\dot{\psi}(\mathbf{r},t) = \int c^2 \cdot \nabla^2 \psi \ dt
$$

we get the momentum of the field and finally:

$$
\psi(\mathbf{r},t) = \int  \dot{\psi}(\mathbf{r},t) dt = \int \int  c^2 \cdot \nabla^2 \psi \ dt^2
$$

to get the (electric) wave field.

### Discretization

The wave field and momentum field are discretized into 2-dimensional grids (formally notated as $\Psi_i[\mathbf{n}]: (\mathbb{N}_0^2 \rightarrow \mathbb{F}_b)_i$ and $\dot{\Psi}_i[\mathbf{n}]: (\mathbb{N}_0^2 \rightarrow \mathbb{F}_b)_i$, where $\mathbb{F}_b$ is the floating point set of bit length $b$) that form a series of discrete timesteps $i$ which evolves the field $\Psi_i[\mathbf{n}] \rightarrow \Psi_{i+1}[\mathbf{n}]$. The discrete timesteps have duration $\Delta t$ and spatial distance $h$ (note that the special distance is the same for both dimensions). The initial conditions of $\Psi_0[\mathbf{n}]$ and $\dot{\Psi}_0[\mathbf{n}]$ are known.

Generally, second order partial deriviative can be written as an concatenation of the Laplacian operator:

$$
(\nabla^2 \circ \psi)(\mathbf{r}, t)
$$

The discrete Laplacian operator can be derived using finite differences. We will denote the discrete version as $D^2$ here. $D^2$ can be written as a stencil:

$$  D^2_h = \frac{1}{h^2} \begin{bmatrix}
0 & 1 & 0\\
1 & -4 & 1\\
0 & 1 & 1
\end{bmatrix} $$

We can get the discrete temporal second-order derviative by convolving the stencil over the discrete wave field and multiplying it with $c^2$:

$$
\ddot{\Psi}_{i+1}[\mathbf{n}] = c^2 \cdot (D^2_h * \Psi_i)[\mathbf{n}]
$$

using the leapfrog integration, we can compute the wave field at timestep $i+1$:

$$
\Psi_{i+1}[\mathbf{n}] = \Psi_i[\mathbf{n}] + \dot{\Psi}_i[\mathbf{n}] \cdot \Delta t + \frac{1}{2} \ddot{\Psi}_i[\mathbf{n}] \cdot \Delta t^2
$$

the momentum field at timestep $i+1$ is  computed as

$$
\dot{\Psi}_{i+1}[\mathbf{n}]=\dot{\Psi}_i[\mathbf{n}] + \frac{\ddot{\Psi}_i[\mathbf{n}] + \ddot{\Psi}_{i+1}[\mathbf{n}]}{2} \cdot \Delta t
$$


### Wave Simulation Algorithm

Procedure should be simple:
1. Initialize fields p, pt, ptt
2. For steps i in TimeSet:
    1. Compute p_i+1
    2. compute ptt_i+1 in different buffer
    3. compute pt
    4. display result at step i

---

$$
w = k * v
<=> v = w / k
k = 2 pi / lambda
$$

$$
w = 2 pi / lambda * v
w = 2 pi / T
$$

$$
v = (2 pi / T) / (2 pi / lambda)
v = lambda / T
=> 1.0 = lambda / T
=> we can use whatever T
$$





