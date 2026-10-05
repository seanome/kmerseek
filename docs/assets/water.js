/* The background: one full-screen fragment shader, three things at once.
   Caustic light and rays from the surface (underwater), cells drifting upward with a
   nucleus (inside a tissue), and a slowly turning double helix (inside the sequence).
   It follows the light and dark theme, holds still for prefers-reduced-motion, stops
   when the tab is hidden, and leaves the flat --bg behind it when WebGL is missing. */
(function () {
  'use strict';
  var canvas = document.getElementById('water');
  var btn = document.getElementById('pause');
  var gl = canvas.getContext('webgl', { antialias: false, alpha: false, powerPreference: 'low-power' });
  if (!gl) { canvas.remove(); btn.remove(); return; }

  var VS = 'attribute vec2 p;void main(){gl_Position=vec4(p,0.,1.);}';
  var FS = [
    'precision highp float;',
    'uniform vec2 u_res;uniform float u_t;uniform vec2 u_mouse;uniform float u_light;',
    'float hash(vec2 p){return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}',
    'float noise(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.-2.*f);',
    ' return mix(mix(hash(i),hash(i+vec2(1,0)),f.x),mix(hash(i+vec2(0,1)),hash(i+vec2(1,1)),f.x),f.y);}',
    'float fbm(vec2 p){float v=0.,a=.5;for(int i=0;i<4;i++){v+=a*noise(p);p=p*2.03+vec2(1.7,9.2);a*=.5;}return v;}',
    // caustics: the iterated-sine pattern, warped by noise so it does not tile
    'float caustic(vec2 uv,float t){vec2 w=uv+vec2(fbm(uv*3.+t*.05))*.08;vec2 p=mod(w*6.28318*1.1,6.28318)-250.;vec2 i=p;float c=1.;',
    ' for(int n=0;n<4;n++){float tn=t*(1.-3.5/float(n+1))*.35;i=p+vec2(cos(tn-i.x)+sin(tn+i.y),sin(tn-i.y)+cos(tn+i.x));',
    '  c+=1./length(vec2(p.x/(sin(i.x+tn)*200.),p.y/(cos(i.y+tn)*200.)));}',
    ' c/=4.;c=1.17-pow(c,1.4);return clamp(pow(abs(c),8.),0.,1.);}',
    // one cell: membrane ring, faint body, a nucleus off centre
    'float cell(vec2 uv,vec2 c,float r,float seed){float d=length(uv-c);',
    ' float wob=1.+.06*sin(seed*7.+u_t*.7+atan(uv.y-c.y,uv.x-c.x)*3.);d/=wob;',
    ' float mem=smoothstep(r*.06,0.,abs(d-r))*.9;float body=smoothstep(r,r*.2,d)*.10;',
    ' vec2 nc=c+r*.28*vec2(cos(seed*5.),sin(seed*3.));float nd=length(uv-nc);float nuc=smoothstep(r*.32,r*.12,nd)*.35;',
    ' return mem+body+nuc;}',
    // the double helix seen edge-on: two sine strands, with rungs between them
    'vec3 helix(vec2 uv,float t){float ang=-.55;mat2 R=mat2(cos(ang),-sin(ang),sin(ang),cos(ang));vec2 q=R*(uv-vec2(.62,.42));',
    ' float s=q.x*9.+t*.35,w=q.y;float amp=.075;float a=amp*sin(s),b=amp*sin(s+3.14159);',
    ' float za=cos(s),zb=cos(s+3.14159);',
    ' float sa=smoothstep(.014,.004,abs(w-a))*(.55+.45*za);float sb=smoothstep(.014,.004,abs(w-b))*(.55+.45*zb);',
    ' float per=fract(s/1.05);float rung=smoothstep(.05,.02,abs(per-.5)*1.05)*step(min(a,b),w)*step(w,max(a,b))*(.35+.3*abs(za));',
    ' float fade=smoothstep(.75,.15,abs(q.x))*smoothstep(.34,.1,abs(w));',
    ' return vec3(sa,sb,rung)*fade;}',
    'void main(){vec2 res=u_res;vec2 uv=gl_FragCoord.xy/res;float asp=res.x/res.y;vec2 p=vec2(uv.x*asp,uv.y);',
    ' float t=u_t;p+=(u_mouse-.5)*.03;',
    ' float depth=smoothstep(0.,1.,uv.y+.1*fbm(p*2.+t*.03));',
    ' float ca=caustic(p*.9+vec2(0.,t*.02),t);',
    ' float ray=0.;for(int i=0;i<3;i++){float fi=float(i);float x=p.x*(3.+fi)+fi*1.7+t*(.05+.02*fi)+uv.y*(.6+.3*fi);ray+=pow(max(0.,sin(x)),6.)*(.4-.1*fi);}',
    ' ray*=smoothstep(.2,1.,uv.y);',
    ' float cells=0.;for(int i=0;i<9;i++){float fi=float(i);float h1=hash(vec2(fi,1.)),h2=hash(vec2(fi,2.)),h3=hash(vec2(fi,3.));',
    '  float r=.05+.11*h3;float sp=.012+.02*h1;vec2 c=vec2(h2*asp+.06*sin(t*.3+fi),fract(h1*3.+t*sp)*1.3-.15);',
    '  cells+=cell(p,c,r,fi)*(.5+.5*h3);}',
    ' vec3 hx=helix(p,t);float strands=hx.x+hx.y;',
    // dark water: the marks are light added to a deep background
    ' vec3 dark=mix(vec3(.012,.07,.10),vec3(.03,.22,.28),depth);',
    ' dark+=vec3(.10,.42,.44)*ca*(.25+.35*uv.y);',
    ' dark+=vec3(.35,.75,.75)*ray*.12;',
    ' dark+=vec3(.42,.86,.82)*cells*.32;',
    ' dark+=vec3(1.,.62,.42)*strands*.55+vec3(.9,.8,.6)*hx.z*.45;',
    ' dark*=1.-.35*pow(length(uv-.5),2.2);',
    // sunlit shallows: the same marks, taken out of a pale background so the text stays readable
    ' vec3 lite=mix(vec3(.52,.74,.79),vec3(.88,.95,.955),depth*depth);',
    ' lite+=vec3(.22,.20,.13)*ca*(.30+.40*(1.-uv.y));',
    ' lite+=vec3(.10,.11,.09)*ray;',
    ' lite-=vec3(.10,.17,.15)*cells*.75;',
    ' lite-=vec3(.00,.32,.40)*strands*.55+vec3(.02,.24,.32)*hx.z*.40;',
    ' vec3 col=mix(dark,lite,u_light);',
    ' col+=(hash(gl_FragCoord.xy+t)-.5)*.014;',
    ' gl_FragColor=vec4(col,1.);}'
  ].join('\n');

  function compile(type, src) {
    var s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) { console.error(gl.getShaderInfoLog(s)); return null; }
    return s;
  }
  var vs = compile(gl.VERTEX_SHADER, VS), fs = compile(gl.FRAGMENT_SHADER, FS);
  if (!vs || !fs) { canvas.remove(); btn.remove(); return; }
  var prog = gl.createProgram(); gl.attachShader(prog, vs); gl.attachShader(prog, fs); gl.linkProgram(prog); gl.useProgram(prog);
  gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
  var loc = gl.getAttribLocation(prog, 'p'); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
  var uRes = gl.getUniformLocation(prog, 'u_res'), uT = gl.getUniformLocation(prog, 'u_t'),
      uMouse = gl.getUniformLocation(prog, 'u_mouse'), uLight = gl.getUniformLocation(prog, 'u_light');

  function resize() {
    var d = Math.min(window.devicePixelRatio || 1, 1.5);
    canvas.width = Math.floor(innerWidth * d); canvas.height = Math.floor(innerHeight * d);
    gl.viewport(0, 0, canvas.width, canvas.height);
    if (!running) draw(lastT);   // keep the still frame in step with the new size
  }
  // 1 in light, 0 in dark: the page's own theme, with the data-theme override the explainers use
  var darkQuery = window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)') : null;
  function lightNow() {
    var forced = document.documentElement.getAttribute('data-theme');
    if (forced === 'dark') return 0; if (forced === 'light') return 1;
    return darkQuery && darkQuery.matches ? 0 : 1;
  }
  var mouse = [0.5, 0.5], target = [0.5, 0.5], reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var running = !reduce, frame = null, t0 = performance.now(), lastT = reduce ? 12 : 0;

  function draw(t) {
    lastT = t;
    gl.uniform2f(uRes, canvas.width, canvas.height);
    gl.uniform1f(uT, t);
    gl.uniform2f(uMouse, mouse[0], mouse[1]);
    gl.uniform1f(uLight, lightNow());
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }
  function render(now) {
    mouse[0] += (target[0] - mouse[0]) * .04; mouse[1] += (target[1] - mouse[1]) * .04;
    draw((now - t0) / 1000);
    frame = running && !document.hidden ? requestAnimationFrame(render) : null;
  }
  function setRunning(on) {
    running = on;
    btn.textContent = on ? 'pause the water' : 'let the water move';
    btn.setAttribute('aria-pressed', String(!on));
    if (on && !frame) { t0 = performance.now() - lastT * 1000; frame = requestAnimationFrame(render); }
    if (!on && frame) { cancelAnimationFrame(frame); frame = null; }
  }

  resize();
  window.addEventListener('resize', resize);
  window.addEventListener('pointermove', function (e) { target = [e.clientX / innerWidth, 1 - e.clientY / innerHeight]; }, { passive: true });
  document.addEventListener('visibilitychange', function () {
    if (document.hidden) { cancelAnimationFrame(frame); frame = null; }
    else if (running && !frame) { t0 = performance.now() - lastT * 1000; frame = requestAnimationFrame(render); }
  });
  if (darkQuery && darkQuery.addEventListener) darkQuery.addEventListener('change', function () { if (!running) draw(lastT); });
  btn.addEventListener('click', function () { setRunning(!running); });
  setRunning(running);
  if (!running) draw(lastT);
})();
