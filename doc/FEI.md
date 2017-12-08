<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

#### TARGET_EQUATION
$$f = cos(3x)\*sin(\pi\*y)$$

####@ DOWN_SIDE_BORDER : (j = 0)
$$\underline{u_{i,j}} = 0$$

####@ UP_SIDE_BORDER : (j = 1)
$$\underline{u_{i,j}} = 0$$

####@ LEFT_SIDE_BORDER : (i = 0)
$$\underline{u_{i,j}} = \frac{sin(\pi*y_{j})}{9+{\pi}^2}$$

####@ RIGHT_SIDE_BORDER : (i = \\(\pi\\))
$$\underline{u_{i,j}} = -\frac{sin(\pi*y_{j})}{9+{\pi}^2}$$

####@ FIVE_POINT_DIFFERENCE_FORMAT_EQUATION
$$\underline{u_{i,j}}-\frac{1}{4}\underline{u_{i-1,j}}-\frac{1}{4}\underline{u_{i+1,j}}-\frac{1}{4}\underline{u_{i,j-1}}-\frac{1}{4}\underline{u_{i,j+1}}=\frac{{h}^2}{4}*f(x_{i},y_{j})$$

<meta http-equiv="refresh" content="30">