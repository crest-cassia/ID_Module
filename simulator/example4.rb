#!/home/osaka/.rvm/rubies/ruby-1.9.3-p551/bin/ruby

require 'json'

def norm_rand(mu=0,sigma=1.0)
  r1,r2=Random.rand(),Random.rand();
  Math::sqrt(-1*Math::log(r1))*Math::cos(2*Math::PI*r2)*sigma+mu
end

# x0=ARGV[0].to_f
# x1=ARGV[1].to_f
# x2=ARGV[2].to_f
# x3=ARGV[3].to_f
# x4=ARGV[4].to_f
# x5=ARGV[5].to_f
# x6=ARGV[6].to_f

parsed = JSON.load(open('./_input.json'))

x0 = parsed["x0"]
x1 = parsed["x1"]
x2 = parsed["x2"]
x3 = parsed["x3"]
x4 = parsed["x4"]
x5 = parsed["x5"]
x6 = parsed["x6"]
x7 = parsed["x7"]
x8 = parsed["x8"]
x9 = parsed["x9"]

x10 = parsed["x10"]
x11 = parsed["x11"]
x12 = parsed["x12"]
x13 = parsed["x13"]
x14 = parsed["x14"]
x15 = parsed["x15"]
x16 = parsed["x16"]
x17 = parsed["x17"]
x18 = parsed["x18"]
x19 = parsed["x19"]

x20 = parsed["x20"]
x21 = parsed["x21"]
x22 = parsed["x22"]
x23 = parsed["x23"]
x24 = parsed["x24"]
x25 = parsed["x25"]
x26 = parsed["x26"]
x27 = parsed["x27"]
x28 = parsed["x28"]
x29 = parsed["x29"]

x30 = parsed["x30"]
x31 = parsed["x31"]
x32 = parsed["x32"]
x33 = parsed["x33"]
x34 = parsed["x34"]
x35 = parsed["x35"]
x36 = parsed["x36"]
x37 = parsed["x37"]
x38 = parsed["x38"]
x39 = parsed["x39"]

x40 = parsed["x40"]
x41 = parsed["x41"]
x42 = parsed["x42"]
x43 = parsed["x43"]
x44 = parsed["x44"]
x45 = parsed["x45"]
x46 = parsed["x46"]
x47 = parsed["x47"]
x48 = parsed["x48"]
x49 = parsed["x49"]

x50 = parsed["x50"]
x51 = parsed["x51"]
x52 = parsed["x52"]
x53 = parsed["x53"]
x54 = parsed["x54"]
x55 = parsed["x55"]
x56 = parsed["x56"]
x57 = parsed["x57"]
x58 = parsed["x58"]
x59 = parsed["x59"]

x60 = parsed["x60"]
x61 = parsed["x61"]
x62 = parsed["x62"]
x63 = parsed["x63"]
x64 = parsed["x64"]
x65 = parsed["x65"]
x66 = parsed["x66"]
x67 = parsed["x67"]
x68 = parsed["x68"]
x69 = parsed["x69"]

x70 = parsed["x70"]
x71 = parsed["x71"]
x72 = parsed["x72"]
x73 = parsed["x73"]
x74 = parsed["x74"]
x75 = parsed["x75"]
x76 = parsed["x76"]
x77 = parsed["x77"]
x78 = parsed["x78"]
x79 = parsed["x79"]

x80 = parsed["x80"]
x81 = parsed["x81"]
x82 = parsed["x82"]
x83 = parsed["x83"]
x84 = parsed["x84"]
x85 = parsed["x85"]
x86 = parsed["x86"]
x87 = parsed["x87"]
x88 = parsed["x88"]
x89 = parsed["x89"]

x90 = parsed["x90"]
x91 = parsed["x91"]
x92 = parsed["x92"]
x93 = parsed["x93"]
x94 = parsed["x94"]
x95 = parsed["x95"]
x96 = parsed["x96"]
x97 = parsed["x97"]
x98 = parsed["x98"]
x99 = parsed["x99"]

x100 = parsed["x100"]
x101 = parsed["x101"]
x102 = parsed["x102"]
x103 = parsed["x103"]
x104 = parsed["x104"]
x105 = parsed["x105"]
x106 = parsed["x106"]
x107 = parsed["x107"]
x108 = parsed["x108"]
x109 = parsed["x109"]

x110 = parsed["x110"]
x111 = parsed["x111"]
x112 = parsed["x112"]
x113 = parsed["x113"]
x114 = parsed["x114"]
x115 = parsed["x115"]
x116 = parsed["x116"]
x117 = parsed["x117"]
x118 = parsed["x118"]
x119 = parsed["x119"]

x120 = parsed["x120"]
x121 = parsed["x121"]
x122 = parsed["x122"]
x123 = parsed["x123"]
x124 = parsed["x124"]
x125 = parsed["x125"]
x126 = parsed["x126"]

result = { "y" => 100*x63*100*x117+100*x115*100*x64+100*x17*100*x29+100*x16*100*x89+100*x47*100*x84*100*x21+100*x92*100*x54+100*x114*100*x22+100*x15*100*x126+100*x62*100*x69+100*x126*100*x26*100*x83+100*x35*100*x12+100*x122*100*x34*100*x62+100*x95*100*x67+100*x124*100*x121*100*x39+100*x50*100*x23+100*x58*100*x44+100*x126*100*x10+100*x45*100*x57+100*x84*100*x12+100*x109*100*x29*100*x33+100*x19*100*x97+100*x11*100*x65*100*x115*100*x31+100*x24*100*x4+100*x45*100*x103+100*x20*100*x101+100*x89*100*x103*100*x111+100*x43*100*x1+100*x36*100*x21+100*x20*100*x60+100*x74*100*x22*100*x63*100*x78+100*x13*100*x21*100*x95+100*x108*100*x3+100*x116*100*x77+100*x20*100*x69+100*x36*100*x11+100*x91*100*x75*100*x105*100*x82+100*x46*100*x54+100*x92*100*x0*100*x65+100*x21*100*x63*100*x98+100*x43*100*x65*100*x66+100*x107*100*x31*100*x115+100*x12*100*x33*100*x72+100*x31*100*x105*100*x22+100*x75*100*x47+100*x2*100*x83*100*x71+100*x93*100*x92+100*x76*100*x78+100*x106*100*x54+100*x119*100*x9+100*x105*100*x8*100*x35+norm_rand(mu=0,sigma=1.0) }

JSON.dump(result, open('./_output.json','w'))
