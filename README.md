# Moldings_cpp

Implementation of the Chen method [1] was implemented. This method segments the mesh into $n$ regions with similar normals.
This is done by creating a graph where the nodes are triangles, which are connected by edges when they are neighbors.
Then, a graph cut energy minimization algorithm is used [2].

[1] - Desai Chen, Pitchaya Sitthi-amorn, Justin T Lan, and Wojciech Matusik. Computing and fabricating multiplanar models. In Computer graphics forum, volume 32, pages 305–315. Wiley Online Library, 2013.

[2] - Yuri Boykov, Olga Veksler, and Ramin Zabih. Fast approximate energy minimization via graph cuts. IEEE Transactions on pattern analysis and machine intelligence, 23(11):1222–1239, 2001.
