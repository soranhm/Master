import matplotlib.pyplot as plt
perm = []
p_dom =[]
v_dom =[]
i_dom = []

p_u = []
v_u = []
i_u = []

v_u2 =[]
for line in open('no_p_9m.txt', 'r'):
    a,b,c,d,e,f,g =  line.split()
    perm.append(float(a))
    p_dom.append(float(b))
    v_dom.append(float(c))
    i_dom.append(float(d))
    p_u.append(float(e))
    v_u.append(float(f))
    i_u.append(float(g))

v_u2.append(0)
for i in range(len(perm)-1):
    v_u2.append(v_u[i+1]-v_u[i])


print(len(perm),len(v_u2))
'''
print(perm)
print('\n')
print(p_dom)
print('\n')
print(v_dom)
print('\n')
print(i_dom)
print('\n')
print(p_u)
print('\n')
print(v_u)
'''
print(max(v_dom)-min(v_dom))


plt.figure()
plt.subplot(311)
plt.plot(perm,p_dom,'o-')
plt.title("Pressure")
plt.xscale('log')
plt.ylabel('porous domain')
plt.subplot(312)
plt.plot(perm,v_dom,'o-')
plt.xscale('log')
plt.ylabel('viscous domain')
plt.subplot(313)
plt.plot(perm,i_dom,'o-')
plt.xscale('log')
plt.xlabel('Permabilty')
plt.ylabel('interface domain')
plt.tight_layout()
#plt.savefig("graf_p.png", dpi=1200)

plt.figure()
plt.subplot(311)
plt.plot(perm,p_u,'o-')
plt.title("Velocity")
plt.xscale('log')
plt.ylabel('porous domain')
plt.subplot(312)
plt.plot(perm,v_u,'o-')
plt.xscale('log')
plt.ylabel('viscous domain')
plt.subplot(313)
plt.plot(perm,i_u,'o-')
plt.xscale('log')
plt.xlabel('Permabilty')
plt.ylabel('interface domain')
#plt.savefig("graf_u.png", dpi=1200)
plt.tight_layout()
plt.show()
'''

x2 = [i for i in range(0,105,5)]
x2.append(115)
y2 = [0,-0.03,0.07,0.04,0.35,0.18,0.09,0.07,0.22,1.3,0.63,0.42,0.28,0.14,0.09,0.25,0.15,0.2,0.09,0.07,0.11,0.13]


def Average(lst):
    return sum(lst) / len(lst)

a = Average(y2)
b =a
print(b/2.)
a = [a for i in range(len(x2))]

#plt.plot(x2,y2,'-.')
#plt.plot(x2,a,label='Average: %.2f'%b)
#plt.ylim(0,10)

#plt.savefig('diameter_size.png',dpi = 1200)




p = perm[:8]
i = i_u[:8]

p2 = perm[8:]
i2 = i_u[8:]


print(p[-1])

plt.subplot(221)
plt.plot(p,i,'-.')
plt.gca().invert_xaxis()
plt.ylabel('Velocity [$\mu m /s$]')
plt.xscale('log')
plt.gca().invert_xaxis()

plt.subplot(222)
plt.plot(p2,i2,'-.')
plt.gca().invert_xaxis()
plt.yticks([])
plt.xscale('log')

plt.subplot(2,2,(3,4))

plt.plot(x2,y2,'-.')
plt.title('Permeability [$\mu m^{2}$]')
plt.ylabel('Velocity [$\mu m /s$]')
plt.xlabel('Time [ s ]')
#plt.plot(perm,v_u,'-.')
#plt.xscale('log')
#plt.savefig("graf_u.png", dpi=1200)
plt.tight_layout(pad=0.2)
plt.savefig("velocity_compare.png", dpi=1200)
'''
plt.plot(p,i,'-.')
#plt.gca().invert_xaxis()
plt.ylabel('Velocity [$\mu m /s$]')
plt.xscale('log')
plt.xlabel('Permeability [$\mu m^{2}$]')
plt.savefig("velocity_devl.png", dpi=1200)
'''
