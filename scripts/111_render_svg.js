const {chromium}=require('playwright');
(async()=>{const b=await chromium.launch();const p=await b.newPage({deviceScaleFactor:3.125});
const svg=require('fs').readFileSync(process.argv[2],'utf8');
await p.setContent(`<html><body style="margin:0">${svg.replace('width="170mm"','width="642px"').replace(/height="[0-9.]+mm"/,'')}</body></html>`);
const el=await p.$('svg'); await el.screenshot({path:process.argv[3]}); 
await p.setContent(`<html><body style="margin:0">${svg}</body></html>`);
await p.pdf({path:process.argv[4], width:'170mm', height:process.argv[5], printBackground:true, pageRanges:'1'});
await b.close();})();
