# Exact check definitions

C is the strict causal matrix. For both examples C=[[0,M],[0,0]] and C^2=0. Write G=MM^T, v for the number of bottoms and m for the number of tops. The program reconstructs M from the cover definitions and computes exact characteristic-polynomial coefficients.

The link-and-interval list comprises all chain counts, all strict interval counts, the Myrheim–Meyer input, ordering fraction, interval probabilities and scaled densities, two reduced actions n-2N0 and n-N0, undirected link degrees and their half-degree histogram, minimal-point heights and their normalized profile, height, and five link-graph characteristic polynomials. These five are D-H, I-D^-1 H, I-2D^-1 H and the latter two restricted to nonisolates. H=C+C^T and D contains its degrees. There are no isolates, so the restricted versions coincide with the full versions.

The seven spectral polynomials are for C+C^T, C-C^T, C^T C-CC^T, the corresponding three for B=-I+2C, and C^2+(C^2)^T. The directed link-degree sequence is also compared. Fixed factors and diagonal shifts relative to the paper's half-scaled operators preserve equality. No perturbations are included.

Let q(u)=det(uI-G). Then the first three polynomials are t^(m-v)q(t^2), (-1)^v t^(m-v)q(-t^2), and (-1)^v t^(m-v)q(t)q(-t). The B polynomials follow by exact affine substitution. The last is t^(m+v).

For the 30-element incidence orders, bottom degree 4 and top degree 2 give link-Laplacian polynomial (t-2)^(m-v)q((t-4)(t-2)). Normalized scale s=1 or 2 gives (t-1)^(m-v)(s^2/8)^v q(8(t-1)^2/s^2).

For the 70-element orders, split M into top-size classes k=1,2,3. Recompute each G_k=a_k I+b_k J, giving (a,b)=(1,0),(5,1),(10,5), with class sizes 7,21,35. Bottom degree is 22. Put D0(t)=product_k(t-k), P_a=(t-22)D0-sum_k a_k D0/(t-k), and P_c=(t-22)D0-sum_k(a_k+7b_k)D0/(t-k). The link polynomial is product_k(t-k)^(m_k-7) P_a^6 P_c. For normalized scale s use (t-1)^56 ((t-1)^2-s^2 a)^6 ((t-1)^2-s^2 c), where a=sum_k a_k/(22k) and c=sum_k(a_k+7b_k)/(22k). All divisions in these displayed formulas cancel exactly. The program uses polynomial products and Fraction arithmetic.

These identities avoid a floating-eigenvalue comparison. Full coefficient lists were also compared with a separate direct full-matrix calculation. The coordinate test uses exact squared differences and checks every ordered pair, including incomparable pairs. The program checks the 30-element subdivision data and the 70-element subset construction. The non-embedding conclusions also require the written circle-order and Radon arguments; no formal proof certificate is asserted.
