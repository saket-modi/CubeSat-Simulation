To-do list:
- RK4 integrator
- Define Molniya orbit properly in README file
- Plots:
    1) Orbit (3D)
    2) Position, velocity, and acceleration for the satellite (2D)
- Performance:
    1) Barnes-Hut algorithm for gravity simulation (nlogn instead of n^2 time complexity)
- Station-keeping logic
- Current program uses a "hacky" code to ensure Tkinter doesn't throw an error while trying to create plots in separate threads (threads apart from main);
- - If two threads are created to display the speed & altitude and the 3D visualization, the backend throws an error;
  - Possible solutions:
  - 1) Run pygame in a separate thread and plot in the main thread
    2) Using multiple cores (experimenting)
