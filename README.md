A better condition to derive the matrix X required to maintain the properties of a normal is to use the cross product formula.

$$
n = t\times b
$$

Then this should be true,

$$
Xn = Mt\times Mb
$$

==//TODO Move this line into a note about matrix multiplication, perhaps with notes on column major vs row major==
```math
$$\begin{bmatrix}
m_{00} & m_{01} & m_{02}\\
m_{10} & m_{11} & m_{12}\\
m_{20} & m_{21} & m_{22}
\end{bmatrix} 
\begin{bmatrix}
t_x \\
t_y \\
t_z 
\end{bmatrix} 
=
\begin{bmatrix}
m_{00} t_x & m_{01} t_y & m_{02} t_z\\
m_{10} t_x & m_{11} t_y & m_{12} t_z\\
m_{20} t_x & m_{21} t_y & m_{22} t_z
\end{bmatrix} 
=
\begin{bmatrix}
m_{00}\\
m_{10}\\
m_{20}
\end{bmatrix} t_x +
\begin{bmatrix}
m_{01}\\
m_{11}\\
m_{21}
\end{bmatrix} t_y +
\begin{bmatrix}
m_{02}\\
m_{12}\\
m_{22}
\end{bmatrix} t_z
$$
```

Let M\[n] be the nth column of the matrix M,

$$
Xn = (t_xM_{[0]}+t_yM_{[1]}+t_zM_{[2]})\times (b_xM_{[0]}+b_yM_{[1]}+b_zM_{[2]})
$$
we can distribute the cross product using the distributive property of cross product,
$$
(u+v)\times w=(u\times w)+(v\times w)
$$
$$
\begin{align}
Xn = &(t_xM_{[0]}\times (b_xM_{[0]}+b_yM_{[1]}+b_zM_{[2]}))+ \\
&(t_yM_{[1]}\times (b_xM_{[0]}+b_yM_{[1]}+b_zM_{[2]}))+\\
&(t_zM_{[2]}\times (b_xM_{[0]}+b_yM_{[1]}+b_zM_{[2]}))
\end{align}
$$
$$
\begin{align}
Xn = &(t_xM_{[0]}\times b_xM_{[0]} + t_xM_{[0]}\times b_yM_{[1]} + t_xM_{[0]}\times b_zM_{[2]})+ \\
&(t_yM_{[1]}\times b_xM_{[0]} + t_yM_{[1]}\times b_yM_{[1]} + t_yM_{[1]}\times b_zM_{[2]})+\\
&(t_zM_{[2]}\times b_xM_{[0]} + t_zM_{[2]}\times b_yM_{[1]} + t_zM_{[2]}\times b_zM_{[2]})
\end{align}
$$
$$
\begin{align}
Xn = &(t_xb_x(M_{[0]}\times M_{[0]}) + t_xb_y(M_{[0]}\times M_{[1]}) + t_xb_z (M_{[0]}\times M_{[2]}))+ \\
&(t_yb_x(M_{[1]}\times M_{[0]}) + t_yb_y(M_{[1]}\times M_{[1]}) + t_yb_z(M_{[1]}\times M_{[2]}))+\\
&(t_zb_x(M_{[2]}\times M_{[0]}) + t_zb_y(M_{[2]}\times M_{[1]}) + t_zb_z(M_{[2]}\times M_{[2]}))
\end{align}
$$
A vector cross product with itself results in $0$ so some of the terms fall out,
$$
\begin{align}
Xn = &t_xb_y(M_{[0]}\times M_{[1]}) + t_xb_z (M_{[0]}\times M_{[2]})+ \\
&t_yb_x(M_{[1]}\times M_{[0]}) + t_yb_z(M_{[1]}\times M_{[2]})+\\
&t_zb_x(M_{[2]}\times M_{[0]}) + t_zb_y(M_{[2]}\times M_{[1]})
\end{align}
$$
Using the anticommutative property of the cross product,
$$
a\times b = -b\times a
$$
$$
\begin{align}
Xn = & t_yb_z(M_{[1]}\times M_{[2]})+ t_zb_y(M_{[2]}\times M_{[1]})+\\
&t_zb_x(M_{[2]}\times M_{[0]}) + t_xb_z (M_{[0]}\times M_{[2]})+ \\
&t_xb_y(M_{[0]}\times M_{[1]}) + t_yb_x(M_{[1]}\times M_{[0]})\\
\end{align}
$$
$$
\begin{align}
Xn = & t_yb_z(M_{[1]}\times M_{[2]})- t_zb_y(M_{[1]}\times M_{[2]})+\\
&t_zb_x(M_{[2]}\times M_{[0]}) - t_xb_z (M_{[2]} \times M_{[0]})+ \\
&t_xb_y(M_{[0]}\times M_{[1]}) - t_yb_x(M_{[0]}\times M_{[1]})\\
\end{align}
$$
Factoring out the common factor results in,
$$
\begin{align}
Xn = &(t_yb_z - t_zb_y)(M_{[1]}\times M_{[2]})+\\
&(t_zb_x - t_xb_z)(M_{[2]}\times M_{[0]})+ \\
&(t_xb_y - t_yb_x)(M_{[0]}\times M_{[1]})\\
\end{align}
$$
We can see the elements of $t \times b$ on the left and a 3x3 matrix $M'$ on the right with each term being one column of $M'$.
$$
M' = 
\begin{bmatrix}
M_{[1]}\times M_{[2]} & M_{[2]}\times M_{[0]} & M_{[0]}\times M_{[1]}
\end{bmatrix} 
$$
$$
M' = 
\begin{bmatrix}
M_{11}M_{22} - M_{21}M_{12} & M_{12}M_{20} - M_{22}M_{10} & M_{10}M_{21}-M_{20}M_{11}\\
M_{21}M_{02} - M_{01}M_{22} & M_{22}M_{00} - M_{02}M_{20} & M_{20}M_{01}-M_{00}M_{21}\\
M_{01}M_{12} - M_{11}M_{02} & M_{02}M_{10} - M_{12}M_{00} & M_{00}M_{11}-M_{10}M_{01}\\
\end{bmatrix} 
$$
$$
M' = 
\begin{bmatrix}
C_{00} & C_{01} & C_{02}\\
C_{10} & C_{11} & C_{12}\\
C_{20} & C_{21} & C_{22}\\
\end{bmatrix} 
$$
$$
Xn = (t\times b) M'
$$
Solving for X, we see that X is equal to the transpose of $M'$.
$$
Xn = (M')^Tn
$$
$$
X=(M')^T
$$
~~$M'$ is known as the cofactor matrix(notated cof(M)) and the transpose of the cofactor matrix $(M')^T$ is known as the adjugate matrix(notated $adj(M)$)~~
$M'$ is known as the adjugate matrix(notated $adj(M)$), which is defined as the transpose of the cofactor matrix(notated $cof(M)$). Transposing $M'$ therefore gives us back the cofactor matrix which is what X equals.
$$
X= cof(M) = adj(M)^T = (cof(M)^T)^T
$$


#### Examining the difference between using the transpose of the inverse and the adjugate transpose.

$$
M^{-1}=\frac{adj(M)}{det(M)}
$$

$$
cof(M)=adj(M)^T
$$

$$
cof(M)=det(M)(M^{-1})^T
$$

This explains why simply using the transpose of the inverse of M was not sufficient since the $det(M)$ contribution is missing. For $det(M) < 0$ the sign would be flipped.  


==//TODO Why does Eric chose this order(12, 20, 01)? Perhaps thinking about how the adjugate  which he is claiming this equals will answer this==
==//TODO Work out how M' equals to the cofactor matrix Dale Weiler is talking about==


----
Eric Lengyel
https://terathon.com/blog/transforming-normals.html

Normals revisited by Dale Weiler(Graphitemaster) 
https://github.com/graphitemaster/normals_revisited?tab=readme-ov-file
