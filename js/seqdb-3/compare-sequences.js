
function main()
{
    console.log("compare_sequences_data", compare_sequences_data);
    show_full_sequences();
}

// --------------------------------------------------------------------------------

function show_full_sequences()
{
    const add_seqence = function(tr, seq) {
        for (let aa of seq) {
            const aa_td = document.createElement("td");
            aa_td.classList.add(`aa${aa}`);
            aa_td.innerHTML = aa;
            tr.appendChild(aa_td)
        }
    };
    
    const tab1 = document.createElement("table");
    for (let group of compare_sequences_data.groups) {
        const tr_group_name = document.createElement("tr");
        tr_group_name.innerHTML = `<td class="group-name" xrowspan="${group.seq.length}">${group.name}</td><td class="seq-id"></td><td class="seq-aas"><table><tr></tr></table></td>`;
        tr_group_name.querySelector(".seq-id").innerHTML = group.seq[0].id;
        add_seqence(tr_group_name.querySelector(".seq-aas tr"), group.seq[0].seq);
        tab1.appendChild(tr_group_name);
        
        // const members = tr1.querySelector(".group-members");
        // for (let seq of group.seq) {
        //     const tr2 = document.createElement("tr");
        //     tr2.innerHTML = `<td class="sequence-id">${seq.id}</td><td><table class="seq-aas"><tr></tr></table></td>`;
        //     members.appendChild(tr2);
        // }
        // tab1.appendChild(tr1);
    }
    document.querySelector("#full-sequences").appendChild(tab1);
}

// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", main);
